#! /usr/bin/env Rscript


library('argparser')

parser = arg_parser(description="Retrieve the N nearest neighbors exomes based on coverage profile")
parser = add_argument(parser, "--reference", help="Reference dataset (BED format)")
parser = add_argument(parser, "--output", default="nearest_neighbors", help="Output filename. List of N first neighbors")

args = parse_args(parser)

if (is.na(args$reference)){
    cat("No reference file provided !\n")
    print(parser)
    q()
}


library('data.table')
library('matrixStats')
library('coop')

fast_cov = function(x){
    1/(NROW(x) -1) * crossprod(scale(x, TRUE , FALSE))
}


data = fread(args$reference, header=T, verbose=F, showProgress=T)

# Remove samples with no variance
cols_tbr = c(1:3, which(colSds(as.matrix(data[, 4:ncol(data)])) == 0) + 3)
data[, (cols_tbr):=NULL]

# Scale data
cols = colnames(data)
data[, (cols) := lapply(.SD, scale), .SDcols=cols]
data = data[, ..cols]

# Compute the PCA
eig = eigen(fast_cov(data))
eigen_values = round(eig$values / sum(eig$values) * 100)

# Use the first 10 eigen values
n_values = min(length(which(eigen_values > 0)), 10)
eigenvalues = head(eigen_values, n_values)
eigenvectors = eig$vectors[, 1:n_values, drop=FALSE]

# Matrix of distances
l = sapply(1:length(eigenvalues), function(i){return(matrix(rep(eigenvectors[, i], each=eigenvalues[i]), ncol=eigenvalues[i], byrow=TRUE))})

df = do.call("cbind.data.frame", l)
d = as.matrix(dist(df))
colnames(d) = rownames(d) = colnames(data)

# Sort samples by distance
neighbors = data.frame(sapply(colnames(d), function(sample){return(names(sort(d[sample,])))}), check.names=F)

# Write dataframe
saveRDS(neighbors, file=paste0(args$output, ".RData"))
