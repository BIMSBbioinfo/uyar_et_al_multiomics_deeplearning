# use ssGSEA approach to score each TCGA sample (using gene expression datasets)
# 
library(singscore)
library(data.table)
source('/data/local/buyar/arcas/subtyping_paper/src/common_functions.R')
geneSetFolder <- '/data/local/buyar/arcas/subtyping_paper/data/GeneSets'
ens2hgnc <- readRDS('/data/local/buyar/arcas/subtyping_paper/data/ens2hgnc.RDS') #id mapping 
GDC_data_dir <- '/data/local/buyar/pancancer_diagnostics/data/GDCdata_prepared' # path to TCGA data 
outDir <- '/data/local/buyar/arcas/subtyping_paper/data/geneset_scores/singscore'

# f: path to gmt format file 
read_msigdb <- function(f) {
  geneSets <- lapply(readLines(f), function(x) {
    unlist(strsplit(x, "\t"))
  })
  #add names
  names(geneSets) <- unlist(lapply(geneSets, function(x) x[1]))
  #remove first two 
  geneSets <- lapply(geneSets, function(x) x[-c(1:2)])
  return(geneSets)
}
# d: path to folder that contains cancerSEA cell state signature files 
read_cancersea <- function(d) {
  files <- dir(d, '.txt$', full.names = T)
  geneSets <- lapply(files, 
                               function(f) {
                                 dt <- data.table::fread(f)
                                 return(unique(dt$EnsembleID))
                               })
  names(geneSets) <- gsub(".txt$", "", basename(files))
  return(geneSets)
}

# f: path to file containing immune gene sets used in xCell package 
read_xCell_geneSets <- function(f, group_by_cell_type = FALSE) {
  l <- readLines(f)
  geneSets <- lapply(l[-1], function(x) {
    res <- unlist(strsplit(gsub(" ", "_", x), "\t"))
    return(setdiff(res[-c(1:2)], ''))
  })
  
  names(geneSets) <- unlist(lapply(l[-1], function(x) 
    gsub(" ", "_", unlist(strsplit(x, "\t"))[1])))
  
  # there are multiple sources for cell types, which can be grouped together 
  if(group_by_cell_type == TRUE) {
    dt <- data.table('name' = names(geneSets), 
                     'celltype' = sapply(strsplit(names(geneSets), "_"), 
                                         function(x) {
                                           paste0(x[1:(length(x)-2)], collapse = '_')
                                         }))
    # for each cell type, get a union of genes that come from multiple sources
    geneSets_by_celltype <- lapply(split(dt, dt$celltype), function(x) {
      unique(unlist(geneSets[x$name]))
    })
    return(geneSets_by_celltype)
  }
  
  return(geneSets)
}

# function to score a gene set 
score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}

message(date(), ' => importing gene sets')
# get msigdb genesets
g1 <- read_msigdb(file.path(geneSetFolder, 'msigdb', 'hallmark_genesets.all.v7.0.symbols.gmt'))
# map hgnc symbols to ensembl ids 
g1 <- map_ens2hgnc(g1, ens2hgnc, "hgnc2ens")

g2 <- read_msigdb(file.path(geneSetFolder, 'msigdb', 'canonical_pathways.v7.1.symbols.gmt'))
g2 <- map_ens2hgnc(g2, ens2hgnc, "hgnc2ens")

# get cancersea gene sets 
g3 <- read_cancersea(file.path(geneSetFolder, 'CancerSEA'))

g4 <- read_xCell_geneSets(file.path(geneSetFolder, 'xCell_cell_signature_genes.tsv'), 
                          group_by_cell_type = TRUE)
# map hgnc symbols to ensembl ids
g4 <- map_ens2hgnc(g4, ens2hgnc, "hgnc2ens")

geneSets <- list("msigdb_hallmarks" = g1, "msigdb_canonical_pathways" = g2, 
                 "cancersea" = g3, "xCell" = g4)

message(date(), ' => reading gex data')

#read GEX tables for each TCGA project 
files <- dir(GDC_data_dir, '^TCGA.*.gex.fpkm.RDS', full.names = T)

gex.list <- pbapply::pblapply(files, function(f) {
  gex <- SummarizedExperiment::assay(readRDS(f))
  # remove normal samples / keep one sample per patient
  gex <- subsetAssayData(gex, unique(sub("(.{12}).+$", "\\1", colnames(gex)))) 
  return(gex)
})
names(gex.list) <- gsub(".gex.fpkm.RDS", "", basename(files))

message(date(), ' => scoring gene sets')

# for each project, get gene set scores 
# rank genes 
scores.list <- pbapply::pbsapply(simplify = F, names(gex.list), function(pr) {
  require(singscore)
  gex <- gex.list[[pr]]
  message(date(), " processing ", pr)
  rankData <- singscore::rankGenes(gex)
  scores <- pbapply::pbsapply(simplify = F, names(geneSets), function(outfile_prefix) {
    require(singscore)
    message(outfile_prefix)
    gs <- geneSets[[outfile_prefix]]
    cl <- parallel::makeCluster(20)
    parallel::clusterExport(cl, varlist = c('score_gene_set', 'gs', 'rankData'), envir = environment())
    scores <- do.call(cbind, pbapply::pblapply(cl = cl, names(gs), function(x) {
      d <- data.frame(score_gene_set(rankData, gs[[x]]), row.names = colnames(rankData))
      colnames(d) <- x
      return(d)
    }))
    parallel::stopCluster(cl)
    return(scores)
  })
  return(scores)
})

message(date(), ' => saving results')

# save scores 
for(x in names(geneSets)) {
  dir.create(file.path(outDir, x))
}

for(pr in names(scores.list)) {
  for(gs in names(scores.list[[pr]])) {
    outfile <- file.path(outDir, gs, paste0(pr, ".scores.csv"))
    cat(outfile, "\n")
    write.csv(scores.list[[pr]][[gs]], file = outfile, quote = F)
  }
}

