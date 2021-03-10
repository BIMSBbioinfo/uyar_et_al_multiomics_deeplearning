library(argparser, quietly = TRUE)
library(data.table)

p <- argparser::arg_parser('Setup assays for the pipeline')
p <- argparser::add_argument(p, "--settings", help="path to settings file for the experiment")
argv <- argparser::parse_args(p)

message(date(), " => importing settings from ",file.path(argv$settings))
settings <- yaml::read_yaml(argv$settings)
# read utility functions
source(settings$utility_script)

assayDir <- settings$assay_output$folder
if(settings$features$use_geneset == TRUE) {
  ens2hgnc <- readRDS(settings$ens2hgnc)
  geneSets <- read_msigdb(settings$geneset)
  geneSets.ensembl <- map_ens2hgnc(geneSets, ens2hgnc, "hgnc2ens")
  GOI <- unique(unlist(geneSets.ensembl))
} else {
  GOI <- NULL
}

surv <- do.call(rbind, readRDS(settings$surv))

# path to folder with prepared TCGA data
dataDir <- settings$GDCdataDir

projects <- settings$projects
projectDict <- sapply(yaml::read_yaml(settings$projectDict), 
                      function(x) gsub(' ', '', unlist(strsplit(x, ','))))

if(!dir.exists(assayDir)) {
  dir.create(assayDir)
}

cl <- parallel::makeCluster(length(projects))
parallel::clusterExport(cl = cl, varlist = c('settings', 'projectDict', 'dataDir', 
                                             'surv', 'assayDir', 'GOI'))
pbapply::pblapply(cl = cl, projects, function(pr) {
  require(data.table)
  source(settings$utility_script)
  message(date(), " => Preparing assays for project ", pr, " with cancer types ", paste(projectDict[[pr]], collapse = ',') , 
          " for data types ",paste(settings$datatypes, collapse = ','))
  dataList <- prepareData(dataDir = dataDir, clin = surv, 
                          correctBatchEffects = settings$correctBatches, 
                          topN = settings$features$topN, GOI = GOI, 
                          TCGA_Disease_IDs = projectDict[[pr]], 
                          datatypes = settings$datatypes, 
                          gex_flavor = settings$datatype_flavors$gex, 
                          cnv_flavor = settings$datatype_flavors$cnv) 
  # save assays 
  lapply(sapply(settings$datatype_combinations, function(x) strsplit(x, ',')), 
         function(datatypes) {
           out <- file.path(assayDir, 
                            paste0(pr,'.',
                                   paste(sort(datatypes), collapse = '_')))
           # text input for maui or python-based tools
           if(settings$assay_output$printCSV == TRUE) {
             write.table(x = do.call(rbind, dataList$assay[datatypes]), 
                         file = paste0(out, '.csv'), quote = F, row.names = T, sep = ',')
           }
           # R object input for mofa and other R-based tools 
           if(settings$assay_output$printRDS == TRUE) {
             saveRDS(dataList$assay[datatypes], paste0(out, '.RDS')) 
           }
         })
  # save batch tables per data type 
  lapply(names(dataList$batch), function(x) {
    out <- file.path(assayDir, paste0(pr, '.', x, '.batch_info.tsv'))
    write.table(x = dataList$batch[[x]], 
                file = out, quote = F, row.names = T, sep = '\t')
  })
})
parallel::stopCluster(cl)


