#! /usr/bin/env Rscript

library(argparser)

# Path to HMZDelFinder R script
# src/HMZDelFinder.R
PATH_HMZDELFINDER <- "/ru-auth/local/home/yseeleuthn/scratch/yseeleuthn/Coverage/hmzdelfinder/HMZDelFinder/src/HMZDelFinder.R"

options(warn=1)

parser = arg_parser(description='Running HMZDelFinder')

parser = add_argument(parser, '--input', help='Input bam directory')
parser = add_argument(parser, '--out', help='Output directory')
parser = add_argument(parser, '--bed', help='BED file')
parser = add_argument(parser, '--threads', help='Number of threads [default: 10]')
parser = add_argument(parser, '--data', help='Use an existing data directory')

args = parse_args(parser)

# Display help if either 'input' or 'out' args are empty
if (is.na(args$input) | is.na(args$out)){
    print(parser)
    q()
}

# Default number of threads
if (is.na(args$threads)){
    args$threads = 10
}

bamdir = paste0(args$input, "/")
outdir = paste0(args$out, "/")

# define working directory:
workDir = getwd()
# set project and data directory 
# replace mainDir with the location you want to store experiment results
mainDir = paste0(workDir, "/", outdir, "/"); if (!file.exists(mainDir)){dir.create(mainDir)}

if (is.na(args$data)){
    dataDir = paste0(mainDir, "data/" , sep=""); if (!file.exists(dataDir)){dir.create(dataDir)} # data directory
} else {
    dataDir = paste0(args$data, "/")
}


# Install missing packages from CRAN
list.of.packages <- c("RCurl", "gdata", "data.table", "parallel", "Hmisc", "matrixStats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# Install missing packages from Bioconductor
# Note: for Windows users: Rsubread has to be installed from file  
biocLitePackages <- c("DNAcopy", "GenomicRanges", "Rsubread") 
new.biocLitePackage <- biocLitePackages[!(biocLitePackages %in% installed.packages()[,"Package"])]
if(length(new.biocLitePackage)) { source("http://bioconductor.org/biocLite.R"); biocLite(new.biocLitePackage)}

# load packages
library(RCurl)
library(data.table)
library(gdata)
library(parallel)
library(Hmisc)
library(matrixStats)
library(DNAcopy)
library(GenomicRanges)
library(Rsubread) 

# load HMZDelFinder source code
# Note: source ("https://....")  does not work on some platforms
source(PATH_HMZDELFINDER)
#eval( expr = parse( text = getURL("https://raw.githubusercontent.com/BCM-Lupskilab/HMZDelFinder/master/src/HMZDelFinder.R") ))

# set/create other paths and identifiers
outputDir <- paste0(mainDir, "out/" , sep=""); if (!file.exists(outputDir)){dir.create(outputDir)} # create output directory
plotsDir <- paste0(mainDir, "plots/" , sep=""); if (!file.exists(plotsDir)){dir.create(plotsDir)} # create output plots directory
aohDir <- paste0(mainDir, "AOH/" , sep=""); if (!file.exists(aohDir)){dir.create(aohDir)} 
aohRDataOut <- paste(mainDir, "AOH/extAOH_small.RData", sep="")	# temprary file to store AOH data


# download BED file
# if this does not work, the file can be downloaded from:
# https://www.dropbox.com/s/1v5jbbm2r809ssy/tgp_hg19.bed.tar.gz?dl=0
# and uncompressed manually into dataDir folder
if (is.na(args$bed)){
    if (!file.exists(paste0(dataDir, "tgp_hg19.bed"))){
        if (file.exists(paste0(dataDir, "tgp_hg19.bed.tar.gz"))){file.remove(paste0(dataDir, "tgp_hg19.bed.tar.gz"))}
        dl_from_dropbox( paste0(dataDir, "tgp_hg19.bed.tar.gz"), "1v5jbbm2r809ssy")
        untar(paste0(dataDir, "tgp_hg19.bed.tar.gz"), exdir = dataDir)
    }
    args$bed <- paste0(dataDir, "tgp_hg19.bed") # Default BED file
}


############
## NOTE 1 ## 
############
# To use own WES data and create RPKM files from BAM files one can use calcRPKMsFromBAMs function.
# e.g:
pathToBams <- paste0(workDir, "/", bamdir)
print(pathToBams)
bamFiles <- paste0(pathToBams, dir(pathToBams, "bam$"))
rpkmDir <- dataDir  # place to store RPKM files
sampleNames <- sapply(strsplit(dir(pathToBams, "bam$"), "[/_]"), function(x){x[1]}) # sample identifiers

calcRPKMsFromBAMs(args$bed, bamFiles, sampleNames, rpkmDir, args$threads)

rpkmFiles = paste0(sampleNames, "_rpkm2.txt")           # list of RPKM file names
rpkmFids <- gsub(".rpkm2.txt", "", rpkmFiles) 			# list of sample identifiers
rpkmPaths <- paste0(paste0(dataDir, ""), rpkmFiles) 	# list of paths to RPKM files


############
## NOTE 2 ## 
############
## In this example, we are not performing AOH filtering and VCF files are not required.
## If one wants to perform AOH filtering than need to prepare two lists:
## vcfPaths - the list of paths to VCF files
## vcfFids - the list of sample identifiers that corresponds to VCF files (same order)
## e.g.:
# vcfFiles <- dir (vcfDir,"vcf.bz2$", recursive=T, include.dirs=FALSE)
# vcfFids <- sapply(strsplit(vcfFiles,"[/\\.]"), function(x){x[2]})
# vcfPaths<- paste0(inputDir,vcfFids,"/", sapply(strsplit(snpFiles,"/"), function(x){x[2]}), sep="")
##


########################################
# THRESHOLDS
#
# See description of HMZDelFinder function for details
########################################
is_cmg <- FALSE 		# only for CMG project - otherwhise use FALSE
lowRPKMthreshold <- 0.65# RPKM threshold  
maxFrequency <- 0.05	    # max frequency of HMZ deletion; default =0.005
minAOHsize <- 1000		# min AOH size
minAOHsig <- 0.45		# min AOH signal threshold
vR_id<-"VR"				# ID from VCF FORMAT indicating the number of variant reads, for other variant callers could be "AD"
tR_id<-"DP"				# ID from VCF FORMAT indicating the number total reads 
filter <- "PASS"		# for other variant callers be  '.'

# running HMZDelFinder
results <- runHMZDelFinder (NULL,		# vcfPaths - paths to VCF files for AOH analysis (not used for 1000 genomes) 
                            NULL,		# vcfFids - sample identifiers corresponding to VCF files  (not used for 1000 genomes) 
                            rpkmPaths, 	# paths to RPKM files 
                            rpkmFids,	# samples identifiers corresponding to RPKM files
                            args$threads,	# number of CPU cores
                            aohRDataOut,# temp file to store AOH data
                            args$bed,	# bed file with target 
                            lowRPKMthreshold, #  RPKM threshold 
                            minAOHsize, # min AOH size
                            minAOHsig,	# min AOH signal threshold
                            is_cmg,		# flag used for CMG specific annotations; TRUE samples are from BHCMG cohort, FALSE otherwhise
                            vR_id, 		# ID for 'the number of variants reads' in VCF FORMAT column (default='VR');
                            tR_id,		# ID for 'the number of total reads' in VCF FORMAT column (default='DP')
                            filter)		# only variants with this value in the VCF FILTER column will be used in AOH analysis 


# saving results in csv files
write.csv(results$filteredCalls, paste0(outputDir,"hmzCalls.csv"), row.names=F )

# plotting deletions
#lapply(1:nrow(results$filteredCalls),function(i){
#           plotDeletion(results$filteredCalls, i, results$bedOrdered, results$rpkmDtOrdered,  lowRPKMthreshold, plotsDir, mainText="")})
