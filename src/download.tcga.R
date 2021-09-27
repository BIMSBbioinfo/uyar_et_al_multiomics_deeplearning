# use TCGAbiolinks package to download open-access data from GDC.
# See vignettes for TCGAbiolinks package:
# http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
# download: http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html

library(TCGAbiolinks)

args <- commandArgs(trailingOnly = TRUE)

#workdir for downloading GDCdata
downloaddir <- args[1]
#output folder for writing prepared data after downloads
outDir <- args[2]

pbapply::pboptions(type = 'txt')

get_Data <- function(project, outDir, dataType, ...) {
  if(dataType == 'mut') {
    outFile <- file.path(outDir, 
                         paste0(project, 
                                '.maf'))
    if(!file.exists(outFile)) {
      maf <- GDCquery_Maf(tumor = gsub('TCGA-', '', project), 
                            pipelines = 'muse', save.csv = FALSE)
      write.table(x = maf, file = outFile, 
                  sep = '\t', 
                  row.names = FALSE, 
                  quote = FALSE)
    }
    return(outFile)
  } else if (dataType == 'clinical') {
    outFile <- file.path(outDir, paste0(project, '.clinical.tsv'))
    if(!file.exists(outFile)) {
      clinical <- GDCquery_clinic(project = project, type = "clinical")
      write.table(x = clinical, file = outFile, sep = '\t', quote = FALSE, row.names = FALSE)
    }
    return(outFile)
  }else {
    outFile <- file.path(outDir, 
                         paste(project, dataType, 'RDS', sep = '.'))
    if(!file.exists(outFile)) {
      query <- GDCquery(project = project, ...)
      GDCdownload(query = query)
      dat <- GDCprepare(query = query)
      saveRDS(object = dat, file = outFile)
    }
    return(outFile)
  } 
}

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]

setwd(downloaddir)

cl <- parallel::makeCluster(2, outfile = 'download.tcga.log')
parallel::clusterExport(cl = cl, varlist = c('get_Data', 'outDir'))

parallel::parLapply(cl = cl, X = projects[1:2], fun = function(pr) {
  require("TCGAbiolinks")
  message(date()," => Processing project:",pr)

  message(date()," => getting gex for project:",pr)
  get_Data(project = pr,
            outDir = outDir,
            dataType = 'gex.fpkm',
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification",
            workflow.type = "HTSeq - FPKM")
   
  message(date()," => getting mut for project:",pr)
  get_Data(project = pr,
            outDir = outDir,
            dataType = 'mut')

  message(date()," => getting clinical data for project:",pr)
  get_Data(project = pr,
           outDir = outDir,
           dataType = 'clinical')
   
  message(date()," => getting meth for project:",pr)
  get_Data(project = pr,
           outDir = outDir,
           dataType = 'meth',
           data.category = "DNA Methylation",
           platform = "Illumina Human Methylation 450")
})
parallel::stopCluster(cl)

