import pandas
import os


bam_list = pandas.read_csv(config['input'], header=None, names=['Path'])
bam_list['Name'] = bam_list['Path'].str.extract('.*/([A-Za-z0-9]*)[-_\.].*$')
bam_list = bam_list.set_index('Name')

# Corresp dictionary between name and path to BAM file
data = bam_list.to_dict('index')

# Directory for per sample coverage
dir_coverage = os.path.join(os.path.dirname(config["output"]), "per_sample_coverage/")



def format_entries(wildcards):
    return " ".join(wildcards)


rule all:
    input:
        config["output"]


rule mosdepth:
    input:
        lambda wildcards: data[wildcards.name]['Path']
    output:
        temp(dir_coverage + "{name}.regions.bed.gz")
    log:
        "logs/mosdepth/{name}.log"
    threads:
        2
    shell:
        "mosdepth --by {config[intervals]} --no-per-base --fast-mode --mapq 10 -t 1 " + dir_coverage + "{wildcards.name} {input} &> {log}"


rule extract_column:
    input:
        dir_coverage + "{name}.regions.bed.gz"
    output:
        dir_coverage + "{name}.col.bed"
    params:
        lambda wildcards: bam_list.loc[wildcards.name]['Path']
    shell:
        """
        echo {params} > {output}
        zcat {input} | awk '{{print $NF}}' >> {output}
        """

rule extract_header:
    input:
        expand(dir_coverage + "{name}.regions.bed.gz", name=list(data)[0])
    output:
        temp("header.txt")
    shell:
        """
        echo -e "chr\tstart\tstop" > {output}
        zcat {input} | cut -f 1-3 >> {output}
        """

rule merge:
    """
    Paste all col.bed files to a single file. 'Paste' command cannot be used directly because 
    of the huge number of files.
    """
    
    input:
        files=expand(dir_coverage + "{name}.col.bed", name=data.keys()),
        header="header.txt"
    output:
        config["output"]
    log:
        "logs/mosdepth/merge_ref_db.log"
    threads:
        8
    run:
        merges = []
        files = []
        for i, input_file in enumerate(input.files):
            files.append(input_file)
            if i > 0 and i % 400 == 0:
                shell("paste " + " ".join(files) + " > merge_" + str(i // 400))
                files = []
                merges.append("merge_" + str(i // 400))
        shell("paste " + ' '.join(files) + " > merge_" + str(i // 400 + 1) + " 2>> {log}")
        merges.append("merge_" + str(i // 400 + 1))

        shell("paste {input.header} " + " ".join(merges) + " > {output}" + " 2>> {log}")
        shell("rm merge_*")
