# prepare data for maui 
# download data from Mariathasan paper
# from http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_0.1.13.tar.gz
args <- commandArgs(trailingOnly = TRUE)

sourceFolder <- args[1] #path to source code of IMvigor210CoreBiologies
geneSetFile <- args[2] # path to file with gene names to be included in the analysis

library(SummarizedExperiment)
# load Rdata objects from the downloaded package
load(file.path(sourceFolder, 'data', 'cds.RData'))

counts <- counts(cds) # raw count table
colData <- pData(cds) # sample metadata 
entrez2hdnc <- data.table(fData(cds)) # gene name mapping from entrez to HGNC symbols

normCounts <- sapply(1:nrow(colData), function(i) {
  counts[,i]/colData$sizeFactor[i]
})
colnames(normCounts) <- colnames(counts)

# filter for most variable hallmark + xCell gene sets:
genes <- unique(readLines(geneSetFile)) 
rownames(normCounts) <- entrez2hdnc[match(rownames(normCounts), entrez_id)]$Symbol
normCounts <- normCounts[rownames(normCounts) %in% genes,]
which(table(rownames(normCounts)) > 1)
# drop row with non-unique names
normCounts <- normCounts[-match(names(which(table(rownames(normCounts)) > 1)), 
                                rownames(normCounts)),]


# create folder to save processed data
dir.create(file.path(getwd(), 'assays_immther'))

# Save all normalized genes as RDS
saveRDS(scale(log(normCounts+1)), file.path(getwd(), 'assays_immther', 'mariathasan.gex.all_genes.normalised.RDS'))

# pick top 5000 genes by sd
normCounts <- normCounts[names(sort(apply(normCounts, 1, sd), decreasing = T)[1:5000]),]

# scale and save count table for analysis by dimension reducction tools
out <- file.path(getwd(), 'assays_immther', "mariathasan.gex")
write.table(x = scale(log(normCounts+1)), 
            file = paste0(out, '.csv'), quote = F, row.names = T, sep = ',')

saveRDS(list('gex' = scale(log(normCounts+1))), paste0(out, '.RDS')) 
write.table(x = colData, 
            file = file.path(getwd(), 'assays_immther', 'mariathasan.colData.csv'), 
            quote = F, row.names = T, sep = ',')






