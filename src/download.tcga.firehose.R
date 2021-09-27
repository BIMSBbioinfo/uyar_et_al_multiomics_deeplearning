# Use RTCGAtoolbox package to download processed TCGA datasets such as CNV - GISTIC

library(RTCGAToolbox)
library(TCGAbiolinks)

args <- commandArgs(trailingOnly = TRUE)

pbapply::pboptions(type = 'txt')

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]

#download location for firehose data
downloadDir <- args[1] 
preparedDataDir <- args[2] 

# download TCGA CNV data processed using GISTIC served at BROAD
lastDate <- RTCGAToolbox::getFirehoseAnalyzeDates(last = 1)
cl <- parallel::makeCluster(11)
parallel::clusterExport(cl, varlist = c('lastDate', 'projects', 'downloadDir'))
pbapply::pblapply(X = gsub('TCGA-', '', projects), 
                  FUN = function(pr) {
                    message("\n", date(), " => Downloading data for project: ",pr)
                    if(!file.exists(file.path(downloadDir, 
                                             paste(lastDate, pr, 
                                                   'all_data_by_genes.txt', 
                                                   sep = '-')))) {
                      g <- RTCGAToolbox::getFirehoseData(dataset = pr,
                                                         gistic2Date = lastDate,
                                                         GISTIC = TRUE,
                                                         destdir = downloadDir)
                    } else {
                      message("Data for project ",pr," already downloaded!")
                    }
                  })
parallel::stopCluster(cl)

#NOTE 05.12.2018: SKCM and LAML didn't return the tables for 
#all_data_by_genes.txt, only returned clinical data. 

# process each file and save as RDS files 
assayFiles <- dir(downloadDir, 'all_data_by_genes.txt$', full.names = T)
gistic.assays <- pbapply::pbsapply(simplify = FALSE,  
                          X = assayFiles, function(f) {
                            dt <- data.table::fread(f, header = T, sep = '\t')
                            dt <- dt[,-c('Locus ID', 'Cytoband')]
                            M <- as.matrix(dt[,-1])
                            rownames(M) <- dt$`Gene Symbol`
                            return(M)
                          })
names(gistic.assays) <- paste0("TCGA-", 
                               sapply(strsplit(basename(names(gistic.assays)), "-"),
                                      function(x) x[2]))

# save to file
pbapply::pblapply(names(gistic.assays), function(pr) {
  saveRDS(gistic.assays[[pr]], file = file.path(preparedDataDir, 
                                                paste0(pr, ".cnv_gistic.RDS")))
})






