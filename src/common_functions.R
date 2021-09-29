# common functions used in different analyses
compute_PCA <- function(exp, topN = 50) {
  top <- head(names(sort(apply(exp, 1, sd), decreasing = T)), 5000)
  M <- t(exp[top,])
  pca <- stats::prcomp(M, rank. = topN)
  return(pca)
}

# given TCGA identifiers, create a table of batch information
# see https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
getBatchTable <- function(TCGA_ids) {
  dt <- data.table::data.table(do.call(rbind, lapply(strsplit(TCGA_ids, "-"), function(x) t(data.frame(x)))))
  colnames(dt) <- c('project', 'TSS', 'participant', 'sample_vial', 'portion_analyte', 'plate', 'center')
  dt$bcr_patient_barcode <- paste(dt$project, dt$TSS, dt$participant, sep = '-')
  return(dt)
}


# given a data.table of samples and batch covariates, 
# simplify the batch table: 
# 1. to keep only common batches (e.g. combine batches with few samples)
# 2. remove batch covariates that only contain a single factor 
process_batch_table <- function(df, min_samples_per_batch = 10) {
  l <- lapply(colnames(df), function(n) {
    message(n)
    x <- df[,get(n),drop = T]
    if(length(unique(x)) == 1) {
      return(NULL)
    }
    rare_batches <- names(which(table(x) < min_samples_per_batch))
    x <- ifelse(x %in% rare_batches, 'other', x)
    if(length(unique(x)) == 1) {
      return(NULL)
    }
    if(sum(table(x) > min_samples_per_batch) < 2) {
      return(NULL)
    }
    x <- data.frame(x, stringsAsFactors = T)
    colnames(x) <- n
    return(x)
  })
  dt <- data.table::data.table(do.call(cbind, l[!unlist(lapply(l, is.null))]))
  if(ncol(dt) == 0) {
    return(NULL)
  }
  return(dt)
}


# given the survival data and a matrix (patients on rows, features on columns), 
# return a fit model 
# colnames in M must match bcr_patient_barcode in surv,
# using PFI as default outcome
fit_RSF <- function(surv, M, ntree = 1000, 
                    duration_col = 'PFI.time', 
                    event_col = 'PFI') {
  require(data.table)
  df <- as.data.frame(M)
  df <- cbind(df, data.frame('status' = surv[match(rownames(df), bcr_patient_barcode)][,get(event_col)], 
                             'time' = surv[match(rownames(df), bcr_patient_barcode)][,get(duration_col)]))
  df <- df[!is.na(df$time),]
  fit <- randomForestSRC::rfsrc(formula = Surv(time, status) ~ .,
                                data = df, ntree = ntree, importance = T)
  Cindex <- round(1 - mean(fit$err.rate, na.rm = T), 2)
  message('Cindex',Cindex)
  return(fit)
}

# given a matrix of numerical variables and a data.frame of clinical variables 
# including time and status
# fit a random survival forest and compute c index
fit_RSF_clin <- function(clin, M, ntree = 1000) {
  require(data.table)
  dat <- merge(clin, M, by = 'row.names')
  rownames(dat) <- dat$Row.names
  dat$Row.names <- NULL
  # convert chars to factors
  dat <- do.call(cbind, sapply(simplify = F, colnames(dat), function(x) {
    if(!is.numeric(dat[,x])) {
      dat[,x] <- as.factor(dat[,x])
    }
    return(dat[,x,drop = F])
  }))
  # fit with/without latent factors
  fit1 <- randomForestSRC::rfsrc(formula = Surv(time, status) ~ .,
                                data = dat, ntree = ntree)
  Cindex1 <- round(1 - mean(fit1$err.rate, na.rm = T), 2)
  fit2 <- randomForestSRC::rfsrc(formula = Surv(time, status) ~ .,
                                 data = subset(dat, select = colnames(clin)), 
                                 ntree = ntree)
  Cindex2 <- round(1 - mean(fit2$err.rate, na.rm = T), 2)
  return(list('fit_with' = fit1, 'fit_wo' = fit2, 'C_with' = Cindex1, 'C_wo' = Cindex2))
}

# Use RUVseq to correct for batch effects while preserving tumor/normal status
# counts: raw count matrix
# colData: batch table including at least one column that contains the factor of interest to be preserved
# preserved_factor: variable of interest to preserve (e.g. tumor/normal status)
# k: number of unwanted sources of variation to look for
correct_batch_effects <- function(counts, colData, preserved_factor, k = 10) {
  require(RUVSeq)
  require(EDASeq)
  set <- EDASeq::newSeqExpressionSet(counts = counts,
                             phenoData = colData)
  differences <- RUVSeq::makeGroups(colData[,preserved_factor])
  set_s <- RUVSeq::RUVs(set, unique(rownames(set)), k=k, differences) 
  return(normCounts(set_s))
}





# given a data frame of categorical variables, 
# compute all pairwise cramer's v index for categorical variables
compute_cramers_V <- function(df) {
  do.call(rbind, lapply(1:(ncol(df)-1), function(i) {
    do.call(rbind, lapply((i+1):ncol(df), function(j) {
      if(length(unique(df[,i])) < 2 | length(unique(df[,j])) < 2) {
        data.table::data.table('ref' = colnames(df)[i], 'target' = colnames(df)[j], 
                   'cramersV' = NA, 'chisq_pval' = NA)
      } else {
        cv <- cv.test(df[,i], df[,j])
        data.table::data.table('ref' = colnames(df)[i], 'target' = colnames(df)[j], 
                   'cramersV' = cv$cramersV, 'chisq_pval' = cv$chisq_pval)
      }
    }))
  }))
}

# computing cramer's v index. 
# source: https://www.r-bloggers.com/example-8-39-calculating-cramers-v/
cv.test = function(x,y) {
  t <- chisq.test(x, y, correct=FALSE)
  CV = sqrt(t$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(data.frame('cramersV' = as.numeric(CV), 'chisq_pval' = t$p.value))
}

get_tcga_sample_types <- function(barcodes) {
  st <- list('TP' = 'PRIMARY SOLID TUMOR', 
             'TR' = 'RECURRENT SOLID TUMOR', 
             'TB' = 'Primary Blood Derived Cancer-Peripheral Blood', 
             'TRBM' = 'Recurrent Blood Derived Cancer-Bone Marrow', 
             'TAP' = 'Additional-New Primary', 
             'TM' = 'Metastatic',
             'TAM' = 'Additional Metastatic', 
             'THOC' = 'Human Tumor Original Cells', 
             'TBM' = 'Primary Blood Derived Cancer-Bone Marrow', 
             'NB' = 'Blood Derived Normal', 
             'NT' = 'Solid Tissue Normal', 
             'NBC' = 'Buccal Cell Normal', 
             'NEBV' = 'EBV Immortalized Normal', 
             'NBM' = 'Bone Marrow Normal')
  sampleTypes <- do.call(rbind, lapply(names(st), function(x) {
    b <- TCGAbiolinks::TCGAquery_SampleTypes(barcodes, x)
    message(x, ': found ',length(b), ' samples')
    if(length(b) > 0) {
      return(data.table::data.table('barcode' = b, 'sampletype' = x))
    }
  }))
}

subsetAssayData <- function(M, ids, replaceIDs = TRUE) {
  sample_barcodes <- as.character(colnames(M))
  selected <- sample_barcodes[which(sub("^(.{12}).+$", "\\1", sample_barcodes) %in% ids)]
  normal_samples <- TCGAbiolinks::TCGAquery_SampleTypes(selected, 
                                                        typesample = c('NB', 'NBC', 'NEBV',
                                                                       'NT', 'NBM'))
  selected <- setdiff(selected, normal_samples)
  
  # find patients with multiple samples, pick only one sample per patient
  dt <- data.table::data.table('sample' = selected, 'id' = sub("^(.{12}).+$", "\\1", selected))
  selected <- dt[,head(sample, 1),by = id]$V1
  
  M <- M[,which(as.character(colnames(M)) %in% selected)]
  if(replaceIDs == TRUE) {
    colnames(M) <- dt[,head(sample, 1),by = id]$id
  }
  return(M)
}

get_normal_tcga_samples <- function(sample_ids) {
  normal_samples <- TCGAbiolinks::TCGAquery_SampleTypes(sample_ids, 
                                                        typesample = c('NB', 'NBC', 'NEBV',
                                                                       'NT', 'NBM'))
  return(normal_samples)
}

# train a model to measure accuracy of predicting labels from the given matrix 
# M: samples on the row, features on the column
# each label corresponds to a sample on the row
# labels should correspond to rownames of the matrix
predict_cluster_labels <- function(M, labels, partition = 0.7, ...) {
  require(caret)
  require(ranger)
  df <- data.frame(M)
  labels <- as.factor(as.character(labels))
  # for each factor level, do a one-vs-all classification with probabilities
  auc.dt <- do.call(rbind, lapply(levels(labels), function(x) {
    message(x)
    df$y <- as.factor(ifelse(labels == x, 'one', 'rest'))
    set.seed(3456)
    intrain <- caret::createDataPartition(y = df$y, p= partition, list = FALSE)
    training <- df[intrain,]
    testing <- df[-intrain,]
    if(length(unique(training$y)) == 2 && length(unique(testing$y)) == 2) {
      fit <- ranger::ranger(seed = 42, formula = y ~ ., data = training, probability = T, 
                            ...)
      test.pred <- stats::predict(fit, testing)
      train.auc <- as.numeric(pROC::auc(pROC::roc(training$y, fit$predictions[,'one'])))
      test.auc <- as.numeric(pROC::auc(pROC::roc(testing$y, test.pred$predictions[,'one'])))
      return(data.table::data.table('train_auc' = train.auc, 'test_auc' = test.auc, 'label' = x))
    } else {
      return(data.table::data.table('train_auc' = NA, 'test_auc' = NA, 'label' = x))
    }
  }))
  return(auc.dt)
}

# simplified version of predict_cluster_labels_glmnet
predict_categorical <- function(df, ref_class) {
  require(caret)
  require(ranger)
  # for each factor level, do a one-vs-all classification with probabilities
  df$y <- as.factor(ifelse(df$y == ref_class, 'one', 'rest'))
  set.seed(3456)
  intrain <- caret::createDataPartition(y = df$y, p= 0.6, list = FALSE)
  training <- df[intrain,]
  testing <- df[-intrain,]
  # Build the model using the training set
  set.seed(123)
  require(doParallel)
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)
  ctrl <- trainControl("repeatedcv", repeats = 5, number = 5, sampling = 'down')
  model <- train(
    y ~., data = training, method = "glmnet",
    trControl = ctrl,
    tuneLength = 10
  )
  stopCluster(cl)
  return(model)
}

# returnROC: if set to TRUE, will return a ROC object, when FALSE returns AUC stats 
predict_cluster_labels_glmnet <- function(M, labels, returnModel = FALSE, returnROC = FALSE, 
                                          repeats = NULL, sampling = NULL) {
  require(caret)
  require(ranger)
  df <- data.frame(M)
  labels <- as.factor(as.character(labels))
  # for each factor level, do a one-vs-all classification with probabilities
  res <- sapply(simplify = F, levels(labels), function(x) {
    message(x)
    df$y <- as.factor(ifelse(labels == x, 'one', 'rest'))
    df <- df[!is.na(df$y),]
    set.seed(3456)
    intrain <- caret::createDataPartition(y = df$y, p= 0.6, list = FALSE)
    training <- df[intrain,]
    testing <- df[-intrain,]
    if(length(unique(training$y)) == 2 && length(unique(testing$y)) == 2) {
      # Build the model using the training set
      set.seed(123)
      require(doParallel)
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      ctrl <- trainControl("repeatedcv", number = 5)
      if(!is.null(sampling)) {ctrl$sampling <- sampling}
      if(!is.null(repeats)) {ctrl$repeats <- repeats}
      model <- train(
        y ~., data = training, method = "glmnet",
        trControl = ctrl,
        tuneLength = 10
      )
      stopCluster(cl)
      if(returnModel == TRUE) {return(model)}
      # Best tuning parameter
      message("Best model parameters: alpha=", model$bestTune$alpha, 
              " lambda=", model$bestTune$lambda)
      x.test <- model.matrix(y ~., testing)[,-1]
      test.pred <- predict(model, x.test, type = 'prob')
      r <- pROC::roc(testing$y, test.pred[,'one'])
      test.auc <- as.numeric(pROC::ci.auc(pROC::roc(testing$y, test.pred[,'one'])))
      if(returnROC == TRUE) {
        return(r)
      }
      return(data.table::data.table('auc.lower' = test.auc[1], 'auc' = test.auc[2], 
                        'auc.upper' = test.auc[3], 'label' = x))
    } else {
      if(returnROC == TRUE) {
        return(NULL)
      }
      return(data.table::data.table('auc.lower' = NA, 'auc' = NA, 'auc.upper' = NA, 'label' = x))
    }
  })
  if(returnROC == FALSE) {
    return(do.call(rbind, res))
  }
  return(res)
}

predict_labels_glmnet_multiclass <- function(M, labels, partition = 0.6, nodes = 5) {
  df <- data.frame(M)
  df$y <- as.factor(labels)
  set.seed(3456)
  intrain <- caret::createDataPartition(y = df$y, p= partition, list = FALSE)
  training <- df[intrain,]
  testing <- df[-intrain,]
  # Build the model using the training set
  fit.control <- trainControl(method = "cv", number = 5, 
                              summaryFunction = multiClassSummary, classProbs = TRUE)
  set.seed(123)
  require(doParallel)
  cl <- makePSOCKcluster(nodes)
  registerDoParallel(cl)
  model <- train(
    y ~., data = training, method = "glmnet",
    trControl = fit.control, 
    metric = 'logLoss', #for multi-class prediction
    tuneLength = 10
  )
  stopCluster(cl)
  
  # Best tuning parameter
  message("Best model parameters: alpha=", model$bestTune$alpha, 
          " lambda=", model$bestTune$lambda)
  x.test <- model.matrix(y ~., testing)[,-1]
  test.pred <- predict(model, x.test, type = 'prob')
  y_pred <- apply(test.pred, 1, function(x) colnames(test.pred)[which.max(x)]) 
  res <- cbind(data.frame(test.pred), 
               data.frame('obs' = as.factor(testing$y), 
                          'pred' = as.factor(y_pred)))
  return(res)
}


# train a regression model to predict a target variable (y)
# M: samples on the row, features on the column
predict_numerical <- function(M, y, ...) {
  require(caret)
  require(ranger)
  df <- data.frame(M)
  df$y <- y
  df <- df[!is.na(df$y),]
  intrain <- caret::createDataPartition(y = df$y, p= 0.7, list = FALSE)
  training <- df[intrain,]
  testing <- df[-intrain,]
  fit <- ranger::ranger(seed = 42, formula = y ~ ., data = training, importance = 'permutation', 
                        num.trees = 1000)
  cor.test(training$y, fit$predictions)
  test.pred <- stats::predict(fit, testing)
  cor.test(testing$y, test.pred$predictions)
  
  # get top important variables 
  top <- head(names(sort(ranger::importance(fit), decreasing = T)), 10)
  fit2 <- ranger::ranger(seed = 42, formula = y ~ ., data = training[,c(top, 'y')], importance = 'permutation', 
                         num.trees = 1000)
  test.pred <- stats::predict(fit2, testing)
  cor.test(testing$y, test.pred$predictions)
  
}

evaluate_regression_model <- function(y, y_hat) {
  # Model performance metrics
  data.table::data.table(
    RMSE = RMSE(y_hat, y),
    Rsquare = R2(y_hat, y),
    COR = cor(y_hat, y)
  )
}
  
  
# elastic net model (cross-validate using 10 )
predict_numerical_glmnet <- function(M, y, partition = 0.6) {
  require(glmnet)
  require(caret)
  df <- data.frame(M)
  df$y <- y
  df <- df[!is.na(df$y),]
  set.seed(3456)
  intrain <- caret::createDataPartition(y = df$y, p = partition, list = FALSE)
  training <- df[intrain,]
  testing <- df[-intrain,]
  # Build the model using the training set
  set.seed(123)
  model <- caret::train(
    y ~., data = training, method = "glmnet",
    trControl = trainControl("repeatedcv", repeats = 5, 
                             number = 5, allowParallel = F),
    tuneLength = 10
  )
  # Best tuning parameter
  model$bestTune
  preds_test <- predict(model, model.matrix(y ~., testing)[,-1])
  return(list('model' = model, 'preds' = data.frame('y' = testing$y, 
                                                'y_hat' = preds_test)))
}


# elastic net model (cross-validate using 10 )
predict_labels_glmnet <- function(M, y, returnCoefs = FALSE) {
  require(glmnet)
  require(caret)
  df <- data.frame(M)
  df$y <- as.factor(y)
  df <- df[!is.na(df$y),]
  intrain <- caret::createDataPartition(y = df$y, p = 0.6, list = FALSE)
  training <- df[intrain,]
  testing <- df[-intrain,]
  # Build the model using the training set
  set.seed(123)
  require(doParallel)
  cl <- makePSOCKcluster(10)
  registerDoParallel(cl)
  model <- train(
    y ~., data = training, method = "glmnet",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10
  )
  stopCluster(cl)
  
  # Best tuning parameter
  message("Best model parameters: alpha=", model$bestTune$alpha, 
          " lambda=", model$bestTune$lambda)
  x.test <- model.matrix(y ~., testing)[,-1]
  predictions <- predict(model, x.test)
  stats <- caret::confusionMatrix(predictions, as.factor(testing$y))
  data.table::data.table(t(stats$overall))
}
# using pathway scores per patients and clustering memberships of patients
# calculate classification accuracy for cluster labels
compute_biological_homogeneity <- function(cl, ssgsea) {
    message(date(), " => Computing AUC based on random forest classification model")
    auc_scores <- do.call(rbind, pbapply::pblapply(colnames(cl), function(b) {
      message(b)
      common <- intersect(colnames(ssgsea), rownames(cl))
      auc.ssgsea <- predict_cluster_labels(M = t(ssgsea[, common]), 
                                           labels = cl[common,b,drop = T])
      return(data.frame('method' = b, 
                        'auc.test' = mean(auc.ssgsea$test_auc, na.rm = T), 
                        'auc.train' = mean(auc.ssgsea$train_auc, na.rm = T)))
    }))
}

# sort by coefficient of variation
subset_features_by_coefvar <- function(exp, topN = 100) {
  f <- names(head(sort(apply(exp, 1, function(x) sd(x)/abs(mean(x))), decreasing = T),topN))
  return(exp[f,])
}

subset_features_by_variance <- function(exp, topN = 100) {
  f <- names(head(sort(apply(exp, 1, var), decreasing = T),topN))
  return(exp[f,])
}

get_geneset_membership <- function(geneSets, genes) {
  # initialize matrix
  M <- matrix(rep(0, length(genes) * length(geneSets)), 
              nrow = length(genes), ncol = length(geneSets), 
              dimnames = list(genes, names(geneSets)))
  # for each gene, find which gene sets contain the gene
  for(i in 1:length(geneSets)) {
    M[which(genes %in% geneSets[[i]]), i] <- 1 
  }
  return(M)
}

# given a matrix of measurements (rows: patients, columns: different conditions/covariates)
# and a factor vector (representing cluster membership of each patient), 
# find out if any variables show factor-specific enrichment (factor vs rest)
get_factor_specific_variables <- function(M, factors)  {
  res <- do.call(rbind, pbapply::pblapply(colnames(M), function(x) {
    l <- split(M[,x], factors)
    do.call(rbind, lapply(names(l), function(i) {
      g1 <- l[[i]]
      # skip when too few observations
      if(length(g1) < 5) {
        return(NULL)
      } 
      g2 <- as.vector(unlist(l[setdiff(names(l), i)]))
      t <- wilcox.test(g1, g2, 
                       alternative = 'greater')
      fc <- log2((mean(g1)+0.05) / (mean(g2)+0.05))
      data.frame("variable" = x, "ref_cl" = i, "pval" = t$p.value, "log2fc" = fc)
    }))
  }))
  res$padj <- p.adjust(res$pval, method = 'BH')
  return(data.table::as.data.table(res))
}

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

# geneSets: a list of vectors of ensembl gene ids or hgnc symbols 
# ens2hgnc: data frame of mapping gene ids to gene symbols
# direction: 'ens2hgnc' => maps from gene ids to symbols; 'hgnc2ens' maps symbols to gene ids 
map_ens2hgnc <- function(geneSets, ens2hgnc, direction = 'hgnc2ens') {
  mapped <- NULL
  if(direction == 'hgnc2ens') {
    mapped <- lapply(geneSets, function(x) {
      unique(ens2hgnc[ens2hgnc$hgnc_symbol %in% x,]$ref_gene_id)
    })
  } else if (direction == 'ens2hgnc') {
    mapped <- lapply(geneSets, function(x) {
      unique(ens2hgnc[ens2hgnc$ref_gene_id %in% x,]$hgnc_symbol)
    })
  } else {
    stop("Choose ens2hgnc or hgnc2ens for the direction")
  }
  return(mapped)
}

# given a matrix of measurements (rows: patients, columns: different conditions/covariates)
# and a factor vector (representing cluster membership of each patient), 
# find out which covariates show differential values between pairs of factors 
get_differential_factors <- function(M, factors, Nodes = 10) {
  cl <- parallel::makeForkCluster(nnodes = Nodes)
  parallel::clusterExport(cl = cl, varlist = c('M', 'factors'), envir = environment())
  res <- do.call(rbind, pbapply::pblapply(cl = cl, X = colnames(M), FUN = function(x) {
    l <- split(M[,x], factors)
    do.call(rbind, lapply(names(l), function(i) {
      do.call(rbind, lapply(names(l), function(j) {
        g1 <- l[[i]]
        g2 <- l[[j]]
        t <- wilcox.test(g1, g2, alternative = 'greater')
        fc <- log2((mean(g1)+1) / (mean(g2)+1))
        data.frame("variable" = x, "ref_cl" = i, "target_cl" = j, "pval" = t$p.value, "log2fc" = fc)
      }))
    }))
  }))
  parallel::stopCluster(cl)
  res$padj <- p.adjust(res$pval, method = 'BH')
  # find markers that are really specific for each cluster and differential compared to every other cluster
  res <- data.table::as.data.table(res)
  k <- length(unique(factors))
  # find those markers that are differential compared to all target clusters
  res_sig <- res[padj < 0.05, length(unique(target_cl)), by = c('variable', 'ref_cl')][V1 == (k-1)]
  
  res_sig <- lapply(split(res_sig, res_sig$ref_cl), function(m) {
    m <- merge(m[,1:2], res, by = c('variable', 'ref_cl'))[order(padj)]
    m <- m[ref_cl != target_cl, .SD[which.min(padj)],by = variable][order(padj)][,-3]
  })
  
  return(res_sig)
}

# M: matrix, samples on rows, features on columns 
plot_umap <- function(M, factors, returnData = FALSE) {
  require(umap)
  umap.df <- as.data.frame(umap::umap(M)[['layout']])
  umap.df$group <- factors
  colnames(umap.df)[1:2] <- c('UMAP1', 'UMAP2')
  if(returnData == TRUE) {
    return(umap.df)
  }
  ggplot(umap.df, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(color = group))
}

# get the center coordinate of samples in a matrix group by a factor vector
get_centroid <- function(M, factors) {
  centroids <- data.frame(
      pbapply::pbapply(M, 2, function(x) {
        tapply(X = x, factors, FUN = median)
      })
    )
  colnames(centroids) <- c('x', 'y')
  centroids$text.label <- as.factor(rownames(centroids))
  return(centroids)
}

# M: matrix, samples on rows, features on columns 
plot_tsne <- function(M, factors = NULL, show.labels = F, 
                      jitter = FALSE, returnData = FALSE) {
  require(Rtsne)
  if(jitter == TRUE) {
    M <- jitter(M)
  }
  to_keep <- which(apply(M, 1, sd) != 0)
  if(length(to_keep) < nrow(M)) {
    warning("Removing 0 sd samples N=",nrow(M) - length(to_keep))
    M <- M[to_keep, ]
    factors <- factors[to_keep]
  }
  
  tsne <- Rtsne(M) #as.data.frame(umap::umap(M)[['layout']])
  df <- as.data.frame(tsne$Y)
  rownames(df) <- rownames(M)
  df$group <- factors
  colnames(df)[1:2] <- c('tSNE1', 'tSNE2')
  if(returnData == TRUE) {
    return(df)
  } 
  if(is.null(factors)) {
    return(df)
  }
  p <- ggplot(df, aes(x = tSNE1, y = tSNE2)) + 
    geom_point(aes(color = group))
  if(show.labels == TRUE) {
    centroids <- data.frame(get_basis_matrix(tsne$Y, factors))
    colnames(centroids) <- c('x', 'y')
    centroids$group <- as.factor(rownames(centroids))
    p <- p + geom_text(data = centroids, aes(x = x, y = y, 
                                        label = group), 
                       size = 4)
  }
  return(p)
}

# M: matrix, samples on columns, features on rows
plot_pca <- function(M, factors) {
  pca <- compute_PCA(M, topN = 2)
  df <- data.frame(pca$x, check.names = FALSE, row.names = rownames(pca$x))
  df$group <- factors

  var_exp <- round(diag(cov(pca$x))/sum(diag(cov(pca$x))) * 100, 1)
  ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = group), size = 5, alpha = 0.5) +
    theme_bw(base_size = 8) +
    labs(x = paste0('PC1 (',var_exp[['PC1']],'%)'),
         y = paste0('PC2 (',var_exp[['PC2']],'%)')) 
}

plot_diffusion_map <- function(M, factors, returnData = FALSE) {
  require(diffusionMap)
  dmap <- diffusionMap::diffuse(SNFtool::dist2(M, M)^(1/2))[['X']][,1:2]
  df <- as.data.frame(dmap)
  rownames(df) <- rownames(M)
  df$group <- factors
  colnames(df)[1:2] <- c('DC1', 'DC2')
  if(returnData == TRUE) {
    return(df)
  } 
  ggplot(df, aes(x = DC1, y = DC2)) + 
    geom_point(aes(color = group))
}

# do a kruskall-wallis H-test to find out variables that are enriched for batch factors
find_batch_related_variables <- function(M, batch.df) {
  r <- data.table::data.table(do.call(rbind, lapply(colnames(batch.df), function(batch) {
    message(batch)
    do.call(rbind, pbapply::pblapply(X = colnames(M), function(variable) {
      t <- kruskal.test(M[,variable], batch.df[,batch])
      data.table::data.table('batch' = batch, 'variable' = variable, 'kruskal_wallis_pval' = t$p.value)
    }))
  })))
  return(r)
}

# find features that contribute most to classification of sample clusters
# M: samples on rows, labels: corresponding cluster labels of samples
# topN: Number of top features to return
# ... other parameteres to pass to ranger::ranger function 
get_cluster_contributing_features <- function(M, labels, topN = 100, ...) {
  require(ranger)
  df <- data.frame(M)
  sapply(simplify = F, levels(labels), function(x) {
    message(x)
    df$y <- as.factor(ifelse(labels == x, 'one', 'rest'))
    fit <- ranger::ranger(formula = y ~ ., data = df, importance = 'permutation', ...)
    return(names(head(sort(ranger::importance(fit), decreasing = T), topN)))
  })
}

# kmeans clustering of a given matrix for a list k values 
get_kmeans_subtypes <- function(M, k_vals = 2:11, min_cluster_size = 10, nstart = 1000, iter.max = 30) {
  S <- do.call(cbind, lapply(k_vals, function(k) {
    set.seed(seed = 1234)
    res <- kmeans(M, k, nstart = nstart, iter.max = iter.max)
    return(as.numeric(res$cluster))
  }))
  rownames(S) <- rownames(M)
  colnames(S) <- paste0("kmeans_", k_vals)
  
  #drop clustering schemes that yield clusters smaller than 20
  #we want physiologically relevant clusters
  to_drop <- which(!apply(S, 2, function(x) {
    sum(table(x) < min_cluster_size) == 0
  }))
  if(length(to_drop) == ncol(S)) {
    return(NULL)
  } else {
    return(S[,setdiff(colnames(S), names(to_drop)),drop = F])
  }
}

# compute silhouette score for a matrix with a given partition
compute_silh <- function(M, labels) {
  ss <- summary(cluster::silhouette(x = labels, dist = dist(M)))
  #print(ss)
  ss$avg.width
}

# partition.df: 1 column data.frame with patient barcodes on the row
# first column denotes the group membership (e.g. cluster)
get_median_survival_by_group <- function(surv, partition.df, 
                                         duration_col = 'PFI.time', 
                                         event_col = 'PFI') {
  partition.df <- as.data.frame(cbind(partition.df, 
                                      surv[match(rownames(partition.df), bcr_patient_barcode)][,get(event_col)]))
  partition.df <- cbind(partition.df, surv[match(rownames(partition.df), bcr_patient_barcode)][,get(duration_col)])
  colnames(partition.df) <- c('subtype', event_col, duration_col)
  
  # remove those without survival outcome measured
  partition.df <- partition.df[!is.na(partition.df$PFI.time),]
  sort(tapply(partition.df$PFI.time, as.factor(partition.df$subtype), median),
       decreasing = T)
}

plot_survival <- function(surv, partition.df, 
                          duration_col = 'PFI.time', 
                          event_col = 'PFI', 
                          return_pval = FALSE) {
  partition.df <- as.data.frame(cbind(partition.df, 
                                      surv[match(rownames(partition.df), bcr_patient_barcode)][,get(event_col)]))
  partition.df <- cbind(partition.df, surv[match(rownames(partition.df), bcr_patient_barcode)][,get(duration_col)])
  colnames(partition.df) <- c('subtype', event_col, duration_col)
  
  # remove those without survival outcome measured
  partition.df <- partition.df[!is.na(partition.df[,duration_col]),]
  
  surv_object <-  survival::Surv(time = partition.df[,duration_col], event = partition.df[,event_col])
  fit <- surv_fit(surv_object ~ subtype, partition.df)
  if(return_pval == TRUE) {
    x <- survminer::surv_pvalue(fit, partition.df)
    return(x$pval)
  } else{
    return(ggsurvplot(fit, partition.df, pval = TRUE, risk.table = T, surv.median.line = 'hv')) 
  }
}

evaluate_clusters_ssgsea <- function(cl, ssgsea) {
  message(date(), " => Computing AUC based on random forest classification model")
  auc_scores <- do.call(rbind, pbapply::pblapply(colnames(cl), function(b) {
    message(b)
    common <- intersect(colnames(ssgsea), rownames(cl))
    auc.ssgsea <- predict_cluster_labels(M = t(ssgsea[, common]), 
                                         labels = cl[common,b,drop = T])
    auc.ssgsea$cl_method <- b
    return(auc.ssgsea)
    #return(data.frame('method' = b, 
     #                 'auc.test' = mean(auc.ssgsea$test_auc, na.rm = T)))
  }))
}

# given a matrix of measurements (rows: patients, columns: different conditions/covariates)
# and a factor vector (representing cluster membership of each patient), 
# get a basis matrix that is the mean value per factor of variable
get_basis_matrix <- function(M, factors) {
  B <- pbapply::pbapply(M, 2, function(x) {
    tapply(X = x, factors, FUN = mean)
  })
  return(B)
}

# get a path to a file for TCGA maf file 
# return matrix of # of missense mutations in genes vs sample barcodes
parseMAF <- function(project = 'TCGA-BRCA', 
                     dataDir,
                     var_classes = c('Missense_Mutation', 'Nonsense_Mutation')) {
  f <- file.path(dataDir, paste0(project, '.maf'))
  dat <- fread(f, header = T, sep = '\t')
  if(!is.null(var_classes)) {
    dat <- dat[Variant_Classification %in% var_classes]
  }
  dat$mutationID <- paste(dat$Chromosome, dat$Start_Position, dat$End_Position, dat$Strand, 
                          dat$HGVSc, sep = ':')
  dat <- dat[,length(unique(mutationID)),by = c('Gene', 'Tumor_Sample_Barcode')] #replace Gene with Hugo_Symbol if needed
  
  dat <- data.table::dcast.data.table(dat, Gene ~ Tumor_Sample_Barcode, value.var = 'V1', fill = 0)
  M <- as.matrix(dat[,-1])
  rownames(M) <- dat$Gene
  
  return(M)
}


# dataDir: Folder with TCGA data (output of prepared data) ( = '/data/local/buyar/pancancer_diagnostics/data/GDCdata_prepared/')
# project: Name of TCGA project (e.g TCGA-BRCA)
# topN: top feature count to pick. For mutations, it is by sum of mutations, for others it is by coefvar.
# gex_flavor: can be gex_raw or gex.fpkm
# cnv_flavor: can be "cnv" or "cnv_gistic"
prepareData <- function(dataDir, clin, correctBatchEffects = TRUE, 
                        topN = 1000, GOI = NULL,
                        TCGA_Disease_IDs, datatypes = c('gex', 'cnv', 'mut', 'meth'),
                        gex_flavor, 
                        cnv_flavor, ens2hgnc, cpg2hgnc) {
  patientIDs <- unique(clin[project %in% TCGA_Disease_IDs]$bcr_patient_barcode)
  dataList <- list('assay' = list(), 
                   'batch' = list())
  if('mut' %in% datatypes) {
    # get mutation data
    mut <- lapply(TCGA_Disease_IDs, function(pr) {
      message(date(), " => importing mutation data for ",pr)
      mut <- parseMAF(var_classes = NULL, project = pr, dataDir = dataDir)
      mut <- mut[!is.na(rownames(mut)),]
    })
    # find the union of all genes detected for each project
    all_genes <- unique(unlist(lapply(mut, rownames)))
    
    # for each matrix, find missing genes, add zero rows for missing genes 
    # this is much more efficient than Reduce(merge, by rows)
    mut <- lapply(mut, function(x) {
      missing <- setdiff(all_genes, rownames(x)) 
      message(nrow(x), ' missing ', length(missing))
      m <- matrix(data = rep(0, length(missing) * ncol(x)), nrow = length(missing), ncol = ncol(x))
      rownames(m) <- missing
      rbind(x, m)
    })
    common_features <- Reduce(intersect, lapply(mut, rownames))
    # reorder features and concatenate columns
    mut <- do.call(cbind, lapply(mut, function(x) x[common_features,]))
    if(!is.null(GOI)) {
      # we assume genes have ensemble ids and GOI are gene symbols
      s <- intersect(rownames(mut), ens2hgnc[hgnc_symbol %in% GOI]$ref_gene_id) # keep only GOI
      mut <- mut[s,]
    }
    
    # select genes with most non-zero samples
    selected <- head(names(sort(apply(mut, 1, function(x) sum(x > 0)), decreasing = T)), topN)
    
    # remove normal samples and pick one sample per patient
    mut <- subsetAssayData(mut[selected,], patientIDs, replaceIDs = FALSE)
    rownames(mut) <- paste0('mut.', rownames(mut))
    
    # remove patients that don't have any mutation data
    mut <- mut[,names(which(colSums(mut) > 0))]
    
    # get batch info
    batches <- getBatchTable(colnames(mut))
    batches_processed <- process_batch_table(batches)
    if(!is.null(batches_processed)) {
      batches_processed <- data.frame(batches_processed, row.names = colnames(mut))
    }
    
    # scale by column
    mut <- scale(mut)
    
    # replace sample ids with patient ids
    colnames(mut) <- sub("(.{12}).+", "\\1", colnames(mut))
    dataList$assay <- c(dataList$assay, list('mut' = mut))
    dataList$batch <- c(dataList$batch, list('mut' = batches_processed))
  }
  
  if('gex' %in% datatypes) {
    # get gene expression data
    gex <- lapply(TCGA_Disease_IDs, function(pr) {
      message(date(), " => importing gex data for ",pr)
      gex <- SummarizedExperiment::assay(readRDS(file.path(dataDir, paste0(pr, '.', 
                                                                           gex_flavor, 
                                                                           '.RDS'))))
    })
    common_features <- Reduce(intersect, lapply(gex, rownames))
    gex <- do.call(cbind, lapply(gex, function(x) x[common_features,]))
    
    if (!is.null(GOI)) {
      # we assume genes have ensemble ids and GOI are gene symbols
      s <- intersect(rownames(gex), ens2hgnc[hgnc_symbol %in% GOI]$ref_gene_id) # keep only GOI
      gex <- gex[s,]
    }
    
    gex <- subset_features_by_variance(gex, topN = topN)
    rownames(gex) <- paste0('gex.', rownames(gex))
    
    # get batch table for the samples
    batches <- getBatchTable(colnames(gex))
    batches_processed <- process_batch_table(batches)
    if(!is.null(batches_processed)) {
      batches_processed <- data.frame(batches_processed, row.names = colnames(gex))
      if(correctBatchEffects == TRUE) {
        message("Correcting batch effects for ",paste0(colnames(batches_processed), collapse = ':'))
        # correct batch effects if applicable
        # batch table
        bt <- batches_processed
        # find tumor/normal samples and use it in design matrix for batch correction.
        normal <- get_normal_tcga_samples(rownames(bt))
        bt$status <- ifelse(rownames(bt) %in% normal, 'normal', 'tumor')
        
        gex <- correct_batch_effects(counts = gex, colData = bt, preserved_factor = 'status', k = 10)
      }
    }
    gex <- subsetAssayData(gex, patientIDs)
    gex <- scale(gex)
    dataList$assay <- c(dataList$assay, list('gex' = gex))
    dataList$batch <- c(dataList$batch, list('gex' = batches_processed))
  }
  
  if('cnv' %in% datatypes) {
    # get CNV data
    cnv <- lapply(TCGA_Disease_IDs, function(pr) {
      message(date(), " => importing cnv data for ",pr)
      dat <- readRDS(file.path(dataDir, paste(pr, cnv_flavor, 'RDS', sep = '.')))
      if(cnv_flavor == 'cnv') {
        cnv <- as.matrix(dat[,-c(1:3)])
        rownames(cnv) <-  gsub("\\.[0-9]+$", "", dat$`Gene Symbol`)
      } else if (cnv_flavor == 'cnv_gistic') {
        cnv <- dat
      }
      return(cnv)
    })
    common_features <- Reduce(intersect, lapply(cnv, rownames))
    cnv <- do.call(cbind, lapply(cnv, function(x) x[common_features,]))
    
    if (!is.null(GOI)) {
      # we assume cnv rownames are gene symbols and GOI are also gene symbols
      s <- intersect(rownames(cnv), GOI) # keep only GOI
      cnv <- cnv[s,]
    }
    
    selected <- NULL
    if(cnv_flavor == 'cnv') {
      # get top features # pick features which are most often gained/lost
      selected <- head(names(sort(apply(cnv, 1, function(x) sum(x == 0)))), topN) 
    } else if(cnv_flavor == 'cnv_gistic') {
      # pick top features by variance in gistic scores
      selected <- head(names(sort(apply(cnv, 1, var), decreasing = T)), topN) 
    }
    cnv <- cnv[selected,]
    rownames(cnv) <- paste0('cnv.', rownames(cnv))
    
    # remove normal samples
    cnv <- subsetAssayData(cnv, patientIDs, replaceIDs = FALSE)
    if(cnv_flavor == 'cnv') {
      # remove patients that are all zeros
      cnv <- cnv[,names(which(apply(cnv, 2, function(x) sum(x != 0)) > 0))]
    } else if(cnv_flavor == 'cnv_gistic') {
      # remove patients for whom when data is scaled sd equals to NA
      keep <- names(which(!is.na(apply(scale(cnv), 2, sd))))
      cnv <- cnv[,keep]
    }
    
    # get batch info
    batches <- getBatchTable(colnames(cnv))
    batches_processed <- process_batch_table(batches)
    if(!is.null(batches_processed)) {
      batches_processed <- data.frame(batches_processed, row.names = colnames(cnv))
    }
    cnv <- scale(cnv)
    # replace sample ids with patient ids
    colnames(cnv) <- sub("(.{12}).+", "\\1", colnames(cnv))
    dataList$assay <- c(dataList$assay, list('cnv' = cnv))
    dataList$batch <- c(dataList$batch, list('cnv' = batches_processed))
  }
  if('meth' %in% datatypes) {
    # get CNV data
    meth <- lapply(TCGA_Disease_IDs, function(pr) {
      message(date(), " => importing methylation data for ",pr)
      meth <- SummarizedExperiment::assay(
        readRDS(file.path(dataDir, paste0(pr, '.meth.RDS')))
        )
      # replace NA with 0
      meth[is.na(meth)] <- 0 
      return(meth)
    })
    common_features <- Reduce(intersect, lapply(meth, rownames))
    meth <- do.call(cbind, lapply(meth, function(x) x[common_features,]))
    
    if (!is.null(GOI)) {
      # we assume meth rownames are cpg ids and GOI are gene symbols
      s <- intersect(rownames(meth), cpg2hgnc[hgnc_symbol %in% GOI]$cpg_site) # keep only GOI
      meth <- meth[s,]
    }
    
    # get top features by variance
    meth <- subset_features_by_variance(meth, topN = topN)
    rownames(meth) <- paste0('meth.', rownames(meth))
    
    # remove normal samples
    meth <- subsetAssayData(meth, patientIDs, replaceIDs = FALSE)
    
    # get batch info
    batches <- getBatchTable(colnames(meth))
    batches_processed <- process_batch_table(batches)
    if(!is.null(batches_processed)) {
      batches_processed <- data.frame(batches_processed, row.names = colnames(meth))
    }
    meth <- scale(meth)
    # replace sample ids with patient ids
    colnames(meth) <- sub("(.{12}).+", "\\1", colnames(meth))
    dataList$assay <- c(dataList$assay, list('meth' = meth))
    dataList$batch <- c(dataList$batch, list('meth' = batches_processed))
  }
  
  # keep patients for which there is data for all data types 
  common <- Reduce(intersect, lapply(dataList$assay, colnames))
  
  dataList$assay <- lapply(dataList$assay, function(x) x[,common])
  return(dataList)
}


# run kBET algorithm to estimate batch effects 
run_kBET <- function(M, batch_labels, condition_label = NULL, PCA = FALSE) {
  batch.estimate <- kBET::kBET(M, batch_labels, do.pca = PCA, plot = F)
  if(is.na(batch.estimate)) {
    return(NULL)
  }
  x <- batch.estimate$summary['mean',]
  x$cond <- condition_label
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate',title=paste('batch:', condition_label)) +
    theme_bw() +  
    scale_y_continuous(limits=c(0,1))
  return(list('summary' = x, 'plot' = g))
}


process_clinical_covariates <- function(patientIDs) {
  #collect clinical and survival data and subset by patient IDs - save to a file
  clin <- do.call(rbind, readRDS('/data/local/buyar/pancancer_diagnostics/data/GDCdata_processed2/TCGA.clin.RDS'))
  surv <- do.call(rbind, readRDS('/data/local/buyar/pancancer_diagnostics/data/GDCdata_processed2/TCGA.surv.RDS'))
  # process tumor stage info
  surv$ajcc_pathologic_tumor_stage <- gsub('[ABC]$', '', surv$ajcc_pathologic_tumor_stage)
  # merge clin and surv by patient ids and only keep relevant columns 
  #surv is a cleaned up version of clin, however, clin also has some clinical metadata that is useful. 
  #so surv is the main table, we select some necessary fields from clin. 
  clin <- subset(clin[bcr_patient_barcode %in% patientIDs], 
                 select = c('bcr_patient_barcode', 'cigarettes_per_day', 'weight', 
                                  'alcohol_history', 'alcohol_intensity', 'site_of_resection_or_biopsy',
                                  'bmi', 'years_smoked'))
  surv <- subset(surv[bcr_patient_barcode %in% patientIDs], 
                 select = c('bcr_patient_barcode', 'age_at_initial_pathologic_diagnosis', 
                                  'gender', 'ajcc_pathologic_tumor_stage', 'clinical_stage', 
                                  'histological_type', 'histological_grade', 'menopause_status'))
  clinicalData <- merge(clin, surv, by = 'bcr_patient_barcode')
  
  #add tumor purity data from the paper: https://www.nature.com/articles/ncomms9971
  tumor_purity <- data.table::as.data.table(TCGAbiolinks::Tumor.purity)
  tumor_purity$bcr_patient_barcode <- sub("^(.{12}).+$", "\\1", tumor_purity$Sample.ID)
  #we only use the clinically measured IHC values that are obtained from pathology slides. 
  clinicalData$tumor_purity_IHC <- tumor_purity[match(clinicalData$bcr_patient_barcode, tumor_purity$bcr_patient_barcode)]$IHC
  
  clinicalData$tumor_purity_IHC <- as.numeric(gsub(",", '.', clinicalData$tumor_purity_IHC))
  
  
  # now only keep clinical data that is applicable to the type of cancer 
  possible_covariates <- c('age_at_initial_pathologic_diagnosis', 'cigarettes_per_day', 'weight', 
                           'alcohol_history', 'alcohol_intensity', 
                           'bmi', 'years_smoked', 'site_of_resection_or_biopsy',
                           'gender','ajcc_pathologic_tumor_stage', 'clinical_stage', 
                           'histological_type', 'histological_grade', 'menopause_status',
                           'tumor_purity_IHC')
  
  NA_codes <- c('[Not Available]', '[Not Evaluated]',  '[Not Applicable]', '[Unknown]', '[Discrepancy]')
  
  #select fields that can be used as covariates, meaning at least two class factors with defined values 
  selected_covariates <- paste(do.call(c, sapply(simplify = FALSE, possible_covariates, function(f) {
    x <- clinicalData[,get(f)]
    #replace keywords for undefined to NA
    x[x %in% NA_codes] <- NA

    if(!is.numeric(x)) {
      # find factors with too few observations and collate them into 'other' category
      collated <- names(which(table(x) < 5)) 
      if(length(collated) > 0) {
        message("Collating factors within ", f, " with too few observations into => 'other' category")
        x[x %in% collated] <- 'Other'
      }
      # if there is only one category (excluding 'Other'), drop the factor
      if(length(setdiff(unique(x), c(NA, 'Other'))) == 1){
        message("Discard ",f," single factor")
        return(NULL)
      }
    }
    #if the covariate is NA for more than 20% of patients, discard it
    if(sum(is.na(x)) / length(x) > 0.2) {
      message("Discard ",f," > 20 % NA")
      return(NULL)
    }
    return(f)
  })))
  
  #now subset the clinical data for relevant covariates along with survival info
  clinicalData <- as.data.frame(subset(clinicalData, select = c('bcr_patient_barcode', 
                                                                selected_covariates)))
  rownames(clinicalData) <- clinicalData$bcr_patient_barcode
  clinicalData$bcr_patient_barcode <- NULL
  
  # replace some missing points with NA
  for(col in colnames(clinicalData)) {
    clinicalData[[col]][clinicalData[[col]] %in% NA_codes] <- NA
  }
  
  clinicalData <- do.call(cbind, lapply(colnames(clinicalData), function(f) {
    x <- clinicalData[,f,drop = F]
    if(!is.numeric(x[,1])) {
      collated <- names(which(table(x[,1]) < 5)) 
      clinicalData[x[[1]] %in% collated, f] <- 'Other'
    }
    return(clinicalData[,f,drop =F])
  }))
  
  #factorize categorical variables - 
  #because in python the categorical variables must be encoded with a numerical value
  #categorical <- names(which(sapply(clinicalData, class) == 'character'))
  #for(v in categorical) {
  #  clinicalData[[v]] <- as.numeric(as.factor(clinicalData[[v]])) - 1
  #}
  #replace NA values with 'undefined' 
  #clinicalData <- data.frame(apply(clinicalData, 2, function(x) {
    #x[is.na(x)] <- 'undefined'; return(x)}))
  return(clinicalData)
}


get_tumor_purity_estimates <- function(patientIDs = NULL) {
  dt <- data.table::data.table(TCGAbiolinks::Tumor.purity)
  dt <- data.table('Sample.ID' = dt$Sample.ID, 
             'ESTIMATE' = as.numeric(gsub(",", ".", dt$ESTIMATE)), 
             'ABSOLUTE' = as.numeric(gsub(",", ".", dt$ABSOLUTE)), 
             'LUMP' = as.numeric(gsub(",", ".", dt$LUMP)), 
             'IHC' = as.numeric(gsub(",", ".", dt$IHC)), 
             'CPE' = as.numeric(gsub(",", ".", dt$CPE)), 
             'project' = paste0('TCGA-', dt$Cancer.type), 
             'bcr_patient_barcode' = sub("^(TCGA.{8}).+?$", "\\1", dt$Sample.ID))
  if(is.null(patientIDs)) {
    return(dt)
  }
  return(dt[bcr_patient_barcode %in% patientIDs])
}

# build a cox model for each variable along correcting for existing clinical variables
# return list of covariates that are significant 
get_prognostic_covariates <- function(surv, clin = NULL, explored_covariates, 
                                      event_col = 'PFI', duration_col = 'PFI.time') {
  require(survival)
  res <- do.call(rbind, pbapply::pblapply(colnames(explored_covariates), function(v) {
    #message(v,"\n")
    x <- explored_covariates[,v,drop = F]
    if(!is.null(clin)) {
      dat <- merge(clin, x, by = 'row.names')
      rownames(dat) <- dat$Row.names
      dat$Row.names <- NULL
    } else {
      dat <- data.frame(x)
    }
    dat <- cbind(dat, data.frame('status' = surv[match(rownames(dat), bcr_patient_barcode)][,get(event_col)], 
                               'time' = surv[match(rownames(dat), bcr_patient_barcode)][,get(duration_col)]))
    dat <- dat[!is.na(dat$time),]
    
    cox1 <- coxph(Surv(time, status) ~ ., 
                  data=dat, x = TRUE)
    s <- summary(cox1)
    pval <- as.numeric(s$coefficients[v,][5])
    cindex <- s$concordance[['C']]
    return(data.table::data.table('variable' = v, 'pval' = pval, 'cindex' = cindex))
  }))
  if(!is.null(res)) {
    res <- res[pval < 0.05]
    if(nrow(res) == 0) {
      return(NULL)
    }
    return(res)
  } else {
    return(NULL)
  }
}

get_surv_df <- function(surv, M, event_col = 'PFI', duration_col = 'PFI.time') {
  df <- as.data.frame(M)
  df <- cbind(df, data.frame('status' = surv[match(rownames(df), bcr_patient_barcode)][,get(event_col)], 
                             'time' = surv[match(rownames(df), bcr_patient_barcode)][,get(duration_col)]))
  df <- df[!is.na(df$time),]
  return(df)
}

# apply spectral clustering algorithms available in SNFtool 
# to a given matrix with features on columns, samples on rows
get_spectral_clusters <- function(M, k_vals = 3:7) {
  require(SNFtool)
  d <- SNFtool::dist2(M,M)^(1/2)
  am <- affinityMatrix(d, K = 20, sigma = 0.5)
  df <- data.frame(sapply(k_vals, function(k) {
    spectralClustering(am, k)    
  }), row.names = rownames(am))
  colnames(df) <- paste0('k', k_vals)
  return(df)
}

get_snf_clusters <- function(matrix.list, returnFusedNetworkOnly = FALSE) {
  require(SNFtool)
  common <- Reduce(intersect, lapply(matrix.list, rownames))
  message(date(), ' => getting distance matrix')
  W <- lapply(matrix.list, function(x) {
    d <- SNFtool::dist2(as.matrix(x[common,]),as.matrix(x[common,]))^(1/2)
    am <- affinityMatrix(d, K = 20, sigma = 0.5)
  })
  FN <- SNFtool::SNF(W)
  if(returnFusedNetworkOnly == TRUE){
    return(FN)
  }
  # get best clustering k
  message(date(), ' => looking for best k')
  best <- SNFtool::estimateNumberOfClustersGivenGraph(FN, NUMC = 3:7)
  message(date(), ' => clustering')
  res <- do.call(cbind, lapply(names(best), function(x) {
    group = spectralClustering(FN,best[[x]])    # the final subtypes information
    df <- data.frame('snf' = group, row.names = common)
    colnames(df)[1] <- x
    return(df)
  }))
  return(res)
}

process_tcga_subtypes <- function() {
  subtypes <- data.table::as.data.table(TCGAbiolinks::PanCancerAtlas_subtypes())
  # get patient barcodes from TCGA subtypes table
  subtypes$bcr_patient_barcode <- sub("(^TCGA-..-....).+", "\\1", subtypes$pan.samplesID)
  subtypes$project <- paste0('TCGA-', subtypes$cancer.type)
  subtypes[cancer.type == 'AML']$project <- 'TCGA-LAML'
  subtypes[cancer.type == 'OVCA']$project <- 'TCGA-OV'
  subtypes[cancer.type %in% c('GBM', 'LGG')]$project <- 'glioma'
  subtypes[cancer.type %in% c('COAD', 'STAD', 'ESCA', 'READ')]$project <- 'panGI'
  subtypes[grep('NA', Subtype_Selected)]$Subtype_Selected <- NA
  # drop normal sample subtype information
  normal_samples <- TCGAbiolinks::TCGAquery_SampleTypes(subtypes$pan.samplesID, 
                                                        typesample = c('NB', 'NBC', 'NEBV',
                                                                       'NT', 'NBM'))
  subtypes <- subtypes[!pan.samplesID %in% normal_samples][!is.na(Subtype_Selected)]
  return(subtypes)
}

# df: data.frame with a single numerical column with sample ids in rows 
# function finds the cut point in the first column of df that maximizes survival separation into two
find_survival_cutpoint <- function(df, surv, event_col = 'PFI', duration_col = 'PFI.time') {
    require(maxstat)
    df <- cbind(df, data.frame('status' = surv[match(rownames(df), bcr_patient_barcode)][,get(event_col)], 
                               'time' = surv[match(rownames(df), bcr_patient_barcode)][,get(duration_col)]))
    df <- df[!is.na(df$time),]
    colnames(df)[1] <- 'score'
    st <- maxstat::maxstat.test(Surv(time, status) ~ score, data = df, 
                                smethod="LogRank", pmethod="exactGauss", 
                                abseps=0.01)
    return(as.numeric(st$estimate))
}

# given a Matrix of measurements and partitions of rows as factors
# find out which variables show significant differences across factors 
# one-way anova analysis 
compute_anova <- function(M, factors, return_val = 'pvalue') {
  factors <- as.factor(factors)
  if(length(levels(factors)) > 1) {
    dt <- do.call(rbind, lapply(colnames(M), function(x) {
      val <- M[,x]
      df <- data.frame("value" = as.numeric(val), "labels" = factors)
      rownames(df) <- names(val)
      #linear model 
      mod <- lm(value ~ labels, data = df)
      t <- anova(mod)
      effect_size <- sjstats::eta_sq(t)[['etasq']]
      #stats::TukeyHSD(t)
      return(data.table::data.table("pval" = t$`Pr(>F)`[1], "effect_size"  = effect_size, "variable" = x))
    }))
    dt$padj <- p.adjust(dt$pval, method = 'BH')
    return(dt)
  } else {
    return(NULL)
  }
}

plot_label_specific_vars <- function(M, labels, top_vars_per_label = 1, 
                                     top_vars = NULL, 
                                     padj_threshold = 0.05, 
                                     as_heatmap = FALSE, string_width = 10, 
                                     ...) {
  labels <- droplevels(as.factor(labels)) # drop empty levels
  if(is.null(top_vars)) {
    # get top factors per factor
    dt <- get_factor_specific_variables(M, labels)[order(padj)][padj < padj_threshold]
    dt <- do.call(rbind, lapply(split(dt, dt$ref_cl), 
                                function(x) {head(x, top_vars_per_label)}))
    top_vars <- unique(dt$variable)
  }
  B <- get_basis_matrix(M[,top_vars,drop=F], labels)
  if(as_heatmap == TRUE) {
    p <- pheatmap::pheatmap(t(B), scale = 'row', silent = T, ...)
    return(p)
  }
  df <- melt(B)
  df$Var1 <- stringr::str_wrap(df$Var1, width = string_width)
  ggplot(df, aes(x = Var2, y = Var1)) + 
    geom_point(aes(size = value, color = value, alpha = value)) + 
    scale_color_gradient(low = 'white', high = 'red') + 
    labs(x = '', y = '')
}

# given a matrix, remove variables that are redundant (Based on correlation cut-off)
# samples on the rows, variables on the columns
# downsize_byTopVar : integer value. If set, the top most variable features are first selected.
remove_redundant_variables <- function(M, perc = 99, cutoff = NULL, downsize_byTopVar = NULL) { #, cutoff = 0.7) {
  if(!is.null(downsize_byTopVar)) {
    message(date()," => selecting top ",downsize_byTopVar, " features by variance first")
    M <- t(subset_features_by_variance(t(M), downsize_byTopVar))
  }
  message(date()," => computing pairwise correlations")
  x <- Rfast::cora(M)
  if(is.null(cutoff)) {
    message(date(), " => looking for correlation cutoff...")
    # define cutoff based on correlation distribution, remove those above 99% 
    # sample features for efficiency
    xs <- x
    if(nrow(x) > 5000) {
      s <- sample(1:nrow(x), 5000)
      xs <- x[s,s]
    }
    cutoff <- quantile(xs[upper.tri(xs)], 1:100/100)[[perc]]
    message(date(), " => setting correlation cut-off to ",cutoff, " at ",perc,"th percentile")
  }
  require(caret)
  message(date(), " => subsetting matrix to remove redundant variables")
  hc = findCorrelation(x, cutoff=cutoff, names = F)
  hc = sort(hc)
  M <- M[,-c(hc)]
  return(M)
}

compute_AMI <- function(M, k, target_labels) {
  df <-  data.frame(get_kmeans_subtypes(M, k_vals = k, 0))
  colnames(df) <- 'cluster'
  ami <- aricode::AMI(df$cluster, target_labels)
  return(ami)
}

# geneSetScoresFolder: path to folder with gene set scores
# annotation: subfolder to import
import_geneset_scores <- function(geneSetScoresFolder, annotation = 'msigdb_hallmarks') {
  # use hallmark gene sets to annotate cluster specific functions
  files <- dir(file.path(geneSetScoresFolder, annotation), 
               '.scores.csv', full.names = T)
  scores <- do.call(rbind, sapply(simplify = F, files, function(x) {
    dt <- data.table::fread(x)
    M <- as.matrix(dt[,-1])
    rownames(M) <- dt$V1
    return(M)
  }))
  return(scores)
}

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

# import LFs for tools in inputDir
import_LFs <- function(TOOLS, inputDir) {
  LFs <- sapply(simplify = F, TOOLS, function(folder) {
    message('Importing ',folder)
    files <- dir(file.path(inputDir, folder), pattern = 'factors.csv$', full.names = T)
    LFs <- pbapply::pbsapply(simplify = F, files, function(f) {
      dt <- data.table::fread(f)
      M <- as.matrix(dt[,-1])
      rownames(M) <- dt$V1
      return(M)
    })
    names(LFs) <- unlist(lapply(strsplit(basename(names(LFs)), split = '\\.'), function(x) paste(x[1:2], collapse = '.')))
    return(LFs)
  })
}

# import features weights for tools in inputDir
import_feature_weights <- function(TOOLS, inputDir) {
  FWs <- sapply(simplify = F, TOOLS, function(folder) {
    message('Importing ',folder)
    files <- dir(file.path(inputDir, folder), pattern = 'feature_weights.csv$', full.names = T)
    W <- pbapply::pbsapply(simplify = F, files, function(f) {
      dt <- data.table::fread(f)
      M <- as.matrix(dt[,-1])
      rownames(M) <- gsub("^assay: ", "", dt$V1) # required only for maui
      return(M)
    })
    names(W) <- unlist(lapply(strsplit(basename(names(W)), split = '\\.'), 
                              function(x) paste(x[1:2], collapse = '.')))
    return(W)
  })
  return(FWs)
}

# make a roc plot from a named list of roc objects (pROC::roc)
plot_roc_list <- function(rocs) {
  aucs <- sapply(rocs, pROC::auc)
  p1 <- pROC::ggroc(rocs, size = 1)
  p1$data$name <- paste0(p1$data$name, " (AUC=", round(aucs[p1$data$name], 3),")")
  p1 <- p1 + geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
    theme(legend.title = element_blank(), legend.direction = 'vertical')
  return(p1)
}

# given two single-column data.frames sharing row.names, each column containing character vectors
# plot bi-partite visualisation of labels (useful for comparing how different cluster memberships look)
plot_cluster_comparison <- function(df1, df2) {
  df1 <- data.frame(df1, check.names = F)
  df2 <- data.frame(df2, check.names = F)
  if(sum(is.na(df1[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df1[is.na(df1[,1]),1] <- 'Undefined'
  }
  if(sum(is.na(df2[,1])) > 0) {
    warning("Converting NA values to 'Undefined'")
    df2[is.na(df2[,1]),1] <- 'Undefined'
  }
  df1[,1] <- as.factor(df1[,1])
  df2[,1] <- as.factor(df2[,1])
  
  dt <- data.table(merge(df1[,1,drop=F], df2[,1,drop = F], by = 'row.names'))
  labels <- colnames(dt[,2:3])
  colnames(dt) <- c('rn', 'g1', 'g2')
  dt1 <- dt[order(g1), c('rn', 'g1')]
  dt2 <- dt[order(g2), c('rn', 'g2')]
  dt1$r1 <- 1:nrow(dt1)
  dt2$r2 <- 1:nrow(dt2)
  dt <- merge(dt1, dt2, by = 'rn')[order(rn)]
  ami <- aricode::AMI(dt$g1, dt$g2)
  label_pos_left <- -0.05
  label_pos_right <- 1.05
  # plot segments 
  ggplot(dt[order(g1)]) + 
    geom_point(aes(x = 0, y = r1, color = g1))  +
    geom_point(aes(x = 1, y = r2, color = g2)) + 
    geom_segment(aes(x = 0, xend = 1, y = r1, yend = r2, color = g1), alpha = 0.25) +
    geom_label(data = dt[,median(r1), by = g1], aes(x = label_pos_left, y = V1, label = g1, color = g1), hjust = 1) +
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1], 
                 aes(x = label_pos_left, y = median, xend = 0, yend = max, color = g1)) + 
    geom_segment(data = dt[,list('max' = max(r1), 'min' = min(r1), 'median' = median(r1)), by = g1], 
                 aes(x = label_pos_left, y = median, xend = 0, yend = min, color = g1)) + 
    geom_label(data = dt[,median(r2), by = g2], aes(x = label_pos_right, y = V1, label = g2, color = g2), hjust = 0) +
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2], 
                 aes(x = label_pos_right, y = median, xend = 1, yend = max, color = g2)) + 
    geom_segment(data = dt[,list('max' = max(r2), 'min' = min(r2), 'median' = median(r2)), by = g2], 
                 aes(x = label_pos_right, y = median, xend = 1, yend = min, color = g2)) + 
    annotate('label', x = 0, y = max(dt$r1)+1, label = labels[1], vjust = 0) + 
    annotate('label', x = 1, y = max(dt$r2)+1, label = labels[2], vjust = 0) + 
    ggtitle(label = paste0(labels, collapse = " <-> "), 
            subtitle = paste0("Adjusted Mutual Information: ",round(ami,2))) + 
    theme_minimal() + 
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), 
          legend.position = 'none',
          plot.margin = unit(c(0.1, 1.5, 0, 1.5), units = 'in'),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)) + 
    coord_cartesian(clip = 'off') 
}

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

# given a vector of genes and a list of gene sets, compute which gene sets are over-represented in the given gene pool
# genes: query set of genes
# geneSets: list of gene sets
# background: vector of gene universe (optional. if not provided, union of geneSets is used as background)
compute_geneset_enrichment <- function(genes, geneSets, background = NULL) {
  if(is.null(background)) {
    background <- unique(unlist(geneSets))
  }
  res <- do.call(rbind, lapply(names(geneSets), function(gs) {
    x <- geneSets[[gs]]
    o <- length(intersect(genes, x))
    p <- length(intersect(x, background))/length(background) # probability of overlap
    data.table('term' = gs, 'overlap' = o, 'query_size' = length(genes), 'target_size' = length(x), 
               'expected_overlap' = round(p * length(genes), 2), 
               'prob' = p, 
               'pval' = binom.test(x = o, n = length(genes), p = p, alternative = 'greater')[['p.value']])
  }))
  res$padj <- p.adjust(res$pval, method = 'BH')
  return(res)
}

# function to score a gene set
score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}

# find percent contribution of different datatypes to the top factor specific variables 
# LFs: Latent factors; nested list
# FWs: Feature weights for each LF; nested list 
# tool: pca/maui/mofa
# project: e.g. pancancer
# datatype: e.g. cnv_gex_mut
# factors: factor vector (e.g. cancer types for each sample or subtype info for each sample)
# top_lfs: top LF to pick per factor
# top_featuers: top features to pick per latent factor
plot_datatype_contribution_by_factor <- function(LFs, FWs, tool, project, datatype, 
                                                 factors, top_lfs = 1, top_features = 100) {
  experiment <- paste0(project, '.', datatype)
  L <- LFs[[tool]][[experiment]]
  W <- FWs[[tool]][[experiment]]
  dt <- get_factor_specific_variables(L, factors)
  top_vars <- get_top_by_label(dt[padj < 0.05][order(padj)], 'ref_cl', topN = top_lfs)
  # for each factor, pick top lfs and pick top features per lf and plot composition
  freq <- do.call(rbind, lapply(split(top_vars, top_vars$ref_cl), function(x) {
    l <- unique(x$variable) # one or more lfs 
    # get top features per lf and concatenate
    f <- do.call(c, lapply(l, function(y) {
      head(sort(abs(W[,y]), decreasing = T), top_features)
    }))
    freq <- data.table(table(sub("^(.+?)\\..+$", "\\1", unique(names(f)))))
    freq$factor <- unique(x$ref_cl)
    colnames(freq)[1] <- 'omics'
    return(freq)
  }))
  
  ggplot(freq, aes(x = factor, y = N)) + 
    geom_bar(stat = 'identity', aes(fill = omics), 
             position = 'fill') + coord_flip() + 
    labs(x = '', y = 'Percentage of top features contributing to top latent factors')
}



