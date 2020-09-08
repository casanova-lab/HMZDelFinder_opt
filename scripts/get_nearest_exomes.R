#! /usr/bin/env Rscript

library('argparser')

parser = arg_parser(description="Retrieve the N nearest neighbors exomes based on coverage profile")
parser = add_argument(parser, "--input", help="RDS file storing the dataframe of neighbors")
parser = add_argument(parser, "--output", help="Output filename. List of N first neighbors")
parser = add_argument(parser, "--name", help="Path to the BAM file of the sample")
parser = add_argument(parser, "--number", default=200, help="Number of neighbors to be returned [default: 100]")

args = parse_args(parser)

if (is.na(args$input)){
    print(parser)
    q()
}

if (is.na(args$name)){
    cat("No sample provided!\n")
    print(parser)
    q(status=1)
}

if (args$number <= 0){
    cat("ERROR: Number of closest neighbors (-n|--number) should be a positive integer!\n")
    q(status=1)
}


data = readRDS(args$input)

# Check that the sample of interest is in the database
if (! args$name %in% colnames(data)){
    cat(paste("ERROR: ", args$name, " is not in the database!\n"))
    q(status=1)
}
all_neighbors = data[,args$name]


# Maximum number of neighbors
if (length(all_neighbors) < args$number + 1){
    cat(paste("WARNING: Maximum number of neighbors adjusted to", dim(data)[1] - 1, "\n"))
    args$number = length(all_neighbors) - 1
}

selected_neighbors = head(all_neighbors, args$number+1)
write.table(selected_neighbors, file=args$output, quote=F, col.names=F, row.names=F)
