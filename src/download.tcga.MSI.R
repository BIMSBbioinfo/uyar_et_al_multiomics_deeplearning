# import MSI data for TCGA patients 
# see: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html#Microsatellite_data
library(TCGAbiolinks)

args <- commandArgs(trailingOnly = TRUE)

downloaddir <- args[1]
outDir <- args[2]
setwd(downloaddir)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]


cl <- parallel::makeCluster(20, outfile = 'download.tcga.MSI.log')
parallel::clusterExport(cl = cl, varlist = c('outDir'))
parallel::parLapply(cl = cl, X = projects, fun = function(pr) {
  require("TCGAbiolinks")
  message(date()," => Processing project:",pr)
  query <- tryCatch (
    {
      GDCquery(project = pr, 
                    data.category = "Other",
                    legacy = TRUE,
                    access = "open",
                    data.type = "Auxiliary test")  
    }, error = function(cond) {
      return(NULL)
    }
  )
  if(!is.null(query)) {
    GDCdownload(query)
    msi_results <- GDCprepare_clinic(query, "msi")
    write.table(x = msi_results, quote = F, sep = '\t',
                file = file.path(outDir, paste0(pr, ".MSI.tsv")))
  }
})
parallel::stopCluster(cl)


