# https://www.bioconductor.org/packages/release/bioc/html/omicade4.html
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)
library(omicade4)

Assays <- args[1] # path to RDS object containing a list of assay matrices
outFile <- args[2] # path to output file

A <- readRDS(Assays)

nf <- ifelse(ncol(A[[1]]) > 1000, 50, 20) # round(max(10, sqrt(ncol(A[[1]]))))
message(date(), " => Running MCIA looking for ",nf," factors")
mcoin <- mcia(A, cia.nf=nf)

LF <- mcoin$mcoa$SynVar
rownames(LF) <- gsub("\\.", "-", rownames(LF))

FeatureWeights <- mcoin$mcoa$axis

# save factors
write.table(x = LF, 
            file = outFile,
            quote = F, row.names = T, sep = ',')

# save feature weights
write.table(x = FeatureWeights, 
            file = gsub('.csv$', '.feature_weights.csv', outFile),
            quote = F, row.names = T, sep = ',')

message(date(), " => Finished MCIA")