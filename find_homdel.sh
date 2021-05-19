#! /usr/bin/env bash

set -e
set -x

function usage(){
    prog=$(basename $0)
    echo "
    $prog - Find homozygous deletions for a list of BAM files

    USAGE:
        $prog -i <bam_files.list> -l <intervals.bed>
    OPTIONS:
        -i|--input      *   List of BAM files
        -l|--intervals  *   List of intervals
        -d|--db         *   Path to the database (.RData)
        -o|--output         Path to the output directory [default: .]
        -k|--k              Number of controls [default: 100]
        -t|--threads        Number of threads to use [default: 16]
        --sliding-windows   Use sliding windows to detect deletions shorter than interval size. 
                            Argument value has to be WINDOW_SIZE:OVERLAP (see examples)
        -h|--help           This help

    EXAMPLES:
        $prog -i bams.list -l CCDS.bed
        $prog -i bams.list -l CCDS.bed --sliding_windows 100:50
"
}


OPTIONS=i:d:k:l:o:t:h
LONGOPTS=input:,db:,k:,list:,output:threads:,sliding-windows:,help


! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    exit 2
fi

eval set -- "$PARSED"

DATA_DIR="data"
ref_db="ref_db.RData"
threads=16
k=100
sliding_windows=0
params=""
outptut_dir="."

while true; do
    case "$1" in
        -i|--input)
            input_file="$2"; shift 2 ;;
        -d|--db)
            ref_db="$2"; shift 2 ;;
        -k|--k)
            k="$2"; shift 2 ;;
        -l|--intervals)
            intervals="$2"; shift 2 ;;
        -o|--output)
            output_dir="$2"; shift 2 ;;
        -t|--threads)
            threads="$2"; shift 2 ;;
        --sliding-windows)
            sliding_windows=1; params="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            echo "Programming error"; exit 3 ;;
    esac
done

if [[ -z $input_file ]]; then
    echo "A list of BAM files is required as input !"
    usage
    exit 2
fi


if [[ $sliding_windows == 1 && -z $params ]]; then
    echo "Please specify the size of sliding windows and the overlap "
    exit 2
fi

# If sliding windows, compute the new intervals on the fly
if [[ $sliding_windows == 1 ]]; then
    window_size=$(echo "$params" | cut -f 1 -d ':')  
    overlap=$(echo "$params" | cut -f 2 -d ':')

    $(dirname "${BASH_SOURCE[0]}")/scripts/split_intervals.py --window "$window_size" \
        --overlap "$overlap" "$intervals" > "intervals-$$.bed"

    intervals="intervals-$$.bed"
fi



# Get sample name from its path
function get_sample_name(){
    basename "$1" | sed 's/[-_\.].*//g'
}

function get_nearest_neighbors(){
    # Extract the list of nearest neighbors from the database
    bam_file="$1"
    ref_db="$2"
    k="$3"
    outdir="$4"

    name=$(get_sample_name "$bam_file")

    mkdir -p "$outdir"
    $(dirname "${BASH_SOURCE[0]}")/scripts/get_nearest_exomes.R --input "$ref_db" --name "$bam_file" \
    --number "$k" --output "$outdir/controls.list"
}

function link_bam_files(){
    bam_file="$1"
    name=$(get_sample_name "$bam_file")

    mkdir -p "$output_dir/$name/bams/"
    while read path; do
        ln -fs $(readlink -f $path) "$output_dir/$name/bams/"
    done < "$output_dir/$name/controls.list"
}

while read bam_file; do
    name=$(get_sample_name "$bam_file")

    get_nearest_neighbors "$bam_file" "$ref_db" "$k" "$output_dir/$name"

    link_bam_files "$bam_file"

    # Run HMZDelFinder
    $(dirname "${BASH_SOURCE[0]}")/scripts/run_hmzdelfinder.R --input "$output_dir/$name/bams/" \
    --out "$output_dir/$name" --threads "$threads" --data "$DATA_DIR" --bed "$intervals" || true

    # Get results
    head -1 "$output_dir/$name/out/hmzCalls.csv" > "$output_dir/$name/hmzdelfinder_result.csv" || true
    fgrep -w "$name" "$output_dir/$name/out/hmzCalls.csv" >> "$output_dir/$name/hmzdelfinder_result.csv" || true
done < "$input_file"

# Remove tmp files
rm -f "intervals-$$.bed"
