# HMZDelFinder_opt

HMZDelFinder_opt is a method to detect homozygous and hemizygous (HMZ) deletions, including partial exon deletions in particular, in typical laboratory cohorts of whole exome sequencing (WES) data that are generated over time under different experimental conditions. Specifically, HMZDelFinder_opt selects the K nearest neighbors of the sample of interest based on a PCA-derived coverage distance in order to select, a priori, a reference control set with a coverage profile similar to that of the WES sample studied. In addition, HMZDelFinder_opt uses a sliding window approach to detect HMZ partial exon deletions (deletions spanning less than an exon) that would otherwise be missed.

## Prerequisites

- GNU parallel
- mosdepth (0.2.2 or higher)
- R (v3.4 or higher) 
    - R 'data.table' library
    - R 'matrixStats' library
    - R 'coop' library
- Python 3


## Usage

### Create a database

#### Create a database from scratch using a list of BAM files.

1. Prepare a text file storing the path to your BAM files, one file per line. Absolute or relative
path are accepted.

   Example: 

    **input_bams.list**

    ```
    /path/to/1.bam
    ../2.bam
    ```

2. Create a list of intervals at the BED format. Fourth column is optional.

    **intervals.bed**

    ```
    19      18170418        18170423        IL12RB1
    19      18170704        18170895        IL12RB1
    19      18171932        18172007        IL12RB1
    19      18172991        18173087        IL12RB1
    ```



3. Run the `create_db.sh` script.

    ```bash
    ./create_db.sh -i input_bams.list -l intervals.bed -o ref_db.bed
    ```

    This will create the `ref_db.bed` file storing the coverage profile of all samples and a 
    `ref_db.RData` storing the K nearest neighbors of each sample.


#### Update an existing database with new samples

Samples can be added to an existing database using:

```bash
./create_db.sh -i input_bams.list --update ref_db.bed -o ref_db_new.bed
```

###  Find HMZ deletions, including partial exon deletions

HMZDelFinder_opt relies on HMZDelFinder to detect HMZ deletions.

Please download HMZDelFinder separately at [https://github.com/BCM-Lupskilab/HMZDelFinder](https://github.com/BCM-Lupskilab/HMZDelFinder)
and change the path to HMZDelFinder.R inside `scripts/run_hmzdelfinder.R`


Find homozygous deletions using an interval list:

```bash
find_homdel.sh -i input_bams.list -l intervals.bed --db ref_db.bed
```

Find homozygous deletions using sliding windows on an interval list:

```bash
find_homdel.sh -i input_bams.list -l intervals.bed --sliding-windows 100:50 --db ref_db.bed
```

Here, intervals will be split in windows of 100bp with 50bp overlap to identify HMZ partial exon deletions.

