# https://www.bioconductor.org/packages/release/bioc/html/omicade4.html
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)
library(omicade4)

Assays <- args[1] # path to RDS object containing a list of assay matrices
outFile <- args[2] # path to output file
nFactors <- args[3] # number of factors to compute

A <- readRDS(Assays)
nFactors <- as.numeric(nFactors)
message(date(), " => Running MCIA looking for ",nFactors," factors")
mcoin <- mcia(A, cia.nf=nFactors)

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
