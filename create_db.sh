#! /usr/bin/env bash

function usage(){
    prog=$(basename $0)
    echo "
    $prog - Create a reference database for HMZDelFinder_opt

    USAGE:
        $prog -i <bam_files.list> -l <intervals.bed> -o <output_filename>
    OPTIONS:
        -i|--input      *   List of BAM files to add to the database
        -l|--intervals  *   List of intervals (Not mandatory if using '--update')
        -o|--output         Output filename [default: ref_db.tsv]
        -t|--threads        Number of threads to use [default: 16]
        -u|--update         Update the database using the input list of BAM files.
        -h|--help           This help.

    EXAMPLES:
       $prog -i bams.list -l CCDS.bed -o database.tsv
       $prog -i bams.list --update database.tsv -o new_db.tsv
"
}


# Checks that all filenames in a list have a '.bam' extension
# Returns 1 if true, 0 otherwise
function check_bam_extension(){
    list="$0"

    all_valid=1
    while read filename; do
        if [[ ${filename: -4} != ".bam" ]]; then
            echo "File '$filename' does not have a '.bam' extension!"
            all_valid=0
        fi
    done < "$list"

    echo $all_valid
}




OPTIONS=i:l:o:t:u:h
LONGOPTS=input:,list:,output:,threads:,update:,help


! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    exit 2
fi

eval set -- "$PARSED"

output="ref_db.tsv"
threads=16
update=0

while true; do
    case "$1" in
        -i|--input)
            input_file="$2"; shift 2 ;;
        -l|--intervals)
            intervals="$2"; shift 2 ;;
        -o|--output)
            output="$2"; shift 2 ;;
        -t|--threads)
            threads="$2"; shift 2 ;;
        -u|--update)
            update=1; prev_db="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            echo "Programming error"; exit 3 ;;
    esac
done

if [[ -z "$input_file" ]]; then
    echo "A list of BAM files is required as input !"
    usage
    exit 2
fi

if [[ "$update" == 1 ]]; then
    # Retrieve interval list
    intervals="intervals-$$.bed"
    awk '{OFS="\t"; if (NR>1){print $1, $2, $3}}' "$prev_db" > "$intervals"
fi


if [[ -z $intervals ]]; then
    echo "An interval file is needed"
    exit 2
fi

# Check that all files of the input list have a .bam extension
# Tools do not complain if we provide another filetype, so we have to handle this case very early
if [[ check_bam_extension == 0 ]]; then
    exit 2
fi



# Directory where per sample coverage is stored.
# Can be removed when the database has been generated.
DIR_COVERAGE="per_sample_coverage"
DIR_LOG="log"

mkdir -p "$DIR_LOG" "$DIR_COVERAGE"

# Coverage per sample per interval
# Do not reprocess files already processed
in_db=()
if [[ $update == 1 ]]; then
    in_db=($(head -1 $prev_db | cut -f 4- | tr '\t' '\n'))
fi

while read bam_file; do
    if [[ ! -f $DIR_COVERAGE/$(basename $bam_file .bam).regions.bed.gz ]]; then
        if [[ ! ("${in_db[@]}" =~ "$bam_file") ]]; then
            echo "$bam_file"
        fi
    fi
done < "$input_file" > to_be_processed-$$.list

parallel -j "$threads" mosdepth --by "$intervals" --no-per-base --fast-mode --mapq 10 -t 1 "$DIR_COVERAGE/"{/.} {} :::: "to_be_processed-$$.list"



# Extract coverage column
parallel -j "$threads" "echo {} > $DIR_COVERAGE/{/.}.col && \
             zcat $DIR_COVERAGE/{/.}.regions.bed.gz | \
             awk '{print \$NF}' >> $DIR_COVERAGE/{/.}.col" :::: "$input_file"

# Merge files
# Here we can't use paste on all files directly because of the limit on the number of args in a 
# command line. So we use batches of 500 files.
BATCH_SIZE=500

if [[ "$update" == 1 ]]; then
    cp "$prev_db" "header-$$.txt"
    file="to_be_processed-$$.list"
else
    cat <(echo -e "chr\tstart\tstop") <(cut -f 1-3 "$intervals") > "header-$$.txt"
    file="$input_file"
fi


sed "s/.*\//$DIR_COVERAGE\//g; s/\.[^\.]*$/.col/g" "$file" | \
split -l "$BATCH_SIZE" -d - group-$$-
for group in "group-$$-"*; do 
    paste $(cat $group) > merge-$$-${group##group-$$-}
done
paste "header-$$.txt" merge-$$-* > "$output"

# Remove temporary files
if [[ "$update" == 1 ]]; then
    rm -f "intervals-$$.bed"
fi
rm -rf "header-$$.txt" "group-$$-"* "merge-$$-"* "to_be_processed-$$.list"


# Compute distances between samples
$(dirname "${BASH_SOURCE[0]}")/scripts/compute_nn.R -r "$output" -o "${output%.*}"
