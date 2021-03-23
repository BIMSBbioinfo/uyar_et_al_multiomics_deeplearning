args <- commandArgs(trailingOnly = TRUE)
Assays <- args[1] # path to a CSV file
outFile <- args[2] # path to output file
nFactors <- args[3] # number of components to compute
source('/data/local/buyar/arcas/subtyping_paper/src/common_functions.R')
M <- as.matrix(read.csv(Assays, check.names = FALSE))
PCA <- compute_PCA(M, topN = as.numeric(nFactors))
write.table(x = PCA$x, 
            file = outFile,
            quote = F, row.names = T, sep = ',')
# save feature rotation matrix
write.table(x = PCA$rotation, 
            file = gsub('.csv', '.feature_weights.csv', outFile),
            quote = F, row.names = T, sep = ',')

