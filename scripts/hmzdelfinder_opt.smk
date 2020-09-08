import os
import yaml



INTERVALS = {"default": "/cngpfs/fs1/bioinfo/yseeleuthn/Coverage/simulations/intervals/tgp_hg19.bed",
             "ccds": "/cngpfs/fs1/bioinfo/yseeleuthn/Coverage/CCDSb37_50bp.genes.bed",
             "gencode": "/cngpfs/fs1/bioinfo/yseeleuthn/Coverage/simulations/intervals/Gencode_v29_exons.bed",
             "gencode_sliding_windows100": "/cngpfs/fs1/bioinfo/yseeleuthn/Coverage/simulations/intervals/final.bed"}

config.update({
    "scripts": srcdir(""),
    "intervals_coverage_profile": "/cngpfs/fs1/bioinfo/yseeleuthn/Coverage/ref_db_all/CCDSb37_50bp.genes.bed",
    "intervals": INTERVALS["gencode_sliding_windows100"],
    "data_dir": "data",
    "neighbors": 100
})


BAMS = yaml.full_load(open(config["indivs"]))


rule all:
    input:
        expand("{dataset}/hmz_deletions.csv", dataset=BAMS.keys()),


rule get_nearest_exomes:
    input:
        lambda wildcards: BAMS[wildcards.dataset],
    output:
        "{dataset}/neighbors.list"
    shell:
       "{config[scripts]}/get_nearest_exomes.R --input {config[db]} --name {input} --number {config[neighbors]} --output {output}" 


rule create_bam_directory:
    input:
        "{dataset}/neighbors.list"
    output:
        directory("{dataset}/bams/")
    shell:
        """
        mkdir -p {output}
        while read path; do ln -fs $(readlink -f "$path") {output}; done < {input}
        """


rule hmzdelfinder_selection:
    input:
        directory="{dataset}/bams/",
    output:
        "{dataset}/out/hmzCalls.csv"
    params:
        prefix="{dataset}/"
    threads:
        16
    shell:
        """
        {config[scripts]}/analysis.R --input {input.directory} --out {params.prefix} --threads {threads} --data {config[data_dir]} --bed {config[intervals]}
        """

rule extract:
    input:
        "{dataset}/out/hmzCalls.csv"
    output:
        "{dataset}/hmz_deletions.csv"
    shell:
        """
        head -1 {input} > {output}
        fgrep -w {wildcards.dataset} {input} >> {output} || true
        """
