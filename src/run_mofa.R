library(reticulate)
# run this in interactive mode, turn it off for command-line
# use_python('/home/buyar/.cache/basilisk/1.2.0/MOFA2-1.0.0/mofa_env/bin/python3.7', required = TRUE)
args <- commandArgs(trailingOnly = TRUE)
s <- import("sys")
message("Python version:", s$version)
library(ggplot2)
library(MOFA2)

Assays <- args[1] # path to RDS object containing a list of assay matrices
outFile <- args[2] # path to output file
nFactors <- args[3] # number of factors to compute

A <- readRDS(Assays)
MOFAobject <- create_mofa(A)
modelOptions <- MOFA2::get_default_model_options(MOFAobject)
# number of factors is maximally set to one forth of the total number of samples
modelOptions$num_factors <- min(as.numeric(nFactors), floor(median(sapply(A, ncol))/4)) 
MOFAobject <- prepare_mofa(MOFAobject, model_options = modelOptions)
message(date(), " training MOFA")
MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = TRUE)
message(date(), " => Saving factors to file")
factors <- get_factors(MOFAobject.trained, factors = "all")
write.table(x = factors[[1]], 
            file = outFile,
            quote = F, row.names = T, sep = ',')

# extract weights
weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")
write.table(x = do.call(rbind, weights), 
            file = gsub('.csv$', '.feature_weights.csv', outFile),
            quote = F, row.names = T, sep = ',')

message(date(), ' => Finished mofa!')


