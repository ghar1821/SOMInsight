# Run diffcyt for all time-points in one go ----

library(diffcyt)
library(SummarizedExperiment)
library(data.table)
library(tools)
library(flowCore)
library(Spectre)

library(foreach)
# library(doParallel)
library(doSNOW)

# use modified generate clusters
source("M:/givanna/cmv/code/generate_clusters.R")
source("M:/givanna/cmv/code/run_flowsom_modified.R")

# We divide the data into folds. So we need to identify which fold to use to extract
# features from
fold_files_dir <- "M:/givanna/cmv/loocv_folds/"
setwd(fold_files_dir)
fold_files <- sapply(list.files(), function(x) paste0(fold_files_dir, x))

save_dir <- "M:/givanna/cmv/diffcyt_10x10_20_completeTP_marker_from_LS/"
dir.create(save_dir, recursive = TRUE)

# cores=detectCores()
cl <- makeCluster(5) #not to overload your computer
registerDoSNOW(cl)

foreach(i=1:26, .combine = cbind, 
        .packages = c("diffcyt", "SummarizedExperiment", 
                      "data.table", "tools", "flowCore",
                      "Spectre")) %dopar% {
  
  # use modified generate clusters
  source("M:/givanna/cmv/code/generate_clusters.R")
  source("M:/givanna/cmv/code/run_flowsom_modified.R")
  
  func_markers <- c("HCMV_IE", "HCMV_pp65", "CD25", "CD69", "HLA_DR", 'CD38', "CD279_PD1", "Tim3",
                    "CD278", "Granzyme_B", "CD314_NKG2D", "CD86", "CD274")
  lin_markers <- c("CD57", "CD159c_NKG2C", "CD161", "FOXP3", "CD34", 
                   "CD3", "CD304", "CD163", "CD33", "CD14",
                   "CD127", "CD11C", "CD56", "CD19", "CD45RA", "CD4", "CD8A",
                   "CD16", "TCRgd", "CD27", "CD20", "CD45RO", 
                   "CD197_CCR7", "TCRva7.2", "CDK1")
  other_markers <- c()
  
  timepoints <- c("T1", "T2", "T3", "T4")
  
  # Read data ----
  data_dir <- "M:/givanna/cmv/data_complete_timepoints_only/"
  # dat_per_tp <- rbindlist(dat_per_tp)
  
  
  # We divide the data into folds. So we need to identify which fold to use to extract
  # features from
  # fold_files_dir <- "M:/givanna/cmv/loocv_folds/"
  # setwd(fold_files_dir)
  # fold_files <- sapply(list.files(), function(x) paste0(fold_files_dir, x))
  # 
  save_dir <- "M:/givanna/cmv/diffcyt_10x10_20_completeTP_marker_from_LS/"
  # dir.create(save_dir, recursive = TRUE)
  
  meta_multiplication <- c(0, 20, 40, 60)        
          
  fold_file <- fold_files[i]
  
  fold_split <- fread(fold_file)
  fold_train <- fold_split[fold_split$set == 'train',]
  fold_test <- fold_split[fold_split$set == 'test',]
  
  fold_id <- file_path_sans_ext(basename(fold_file))
  save_dir_fold <- paste0(save_dir, fold_id)
  dir.create(save_dir_fold)
  
  
  message(paste("Processing", fold_id))
  
  for (idx in c(1:4)) {
    tp <- timepoints[idx]
    message(paste("Processing", tp))
    
    # check if this is already done
    setwd(save_dir_fold)
    if (file.exists(paste0("Clustered_trainDat_", tp, ".csv")) &
        file.exists(paste0("Clustered_testDat_", tp, ".csv"))) {
      next
    }
    
    setwd(data_dir)
    dat_for_this_tp <- fread(paste0("CMV_FileName_", tp, "_withMetadata.csv"))
    
    # dat_for_this_tp <- do.subsample(dat_for_this_tp, targets=10000)
    
    dat_train <- dat_for_this_tp[dat_for_this_tp$PatientName %in% fold_train$PatientName,]
    dat_test <- dat_for_this_tp[dat_for_this_tp$PatientName %in% fold_test$PatientName, lin_markers, with=FALSE]
    
    # just for testing
    # dat_train <- do.subsample(dat = dat_train, targets = 10000)
    setwd(save_dir_fold)
    flowsom_res <- run.flowsom.modified(dat = dat_train,
                                        use.cols = lin_markers,
                                        xdim = 10,
                                        ydim = 10,
                                        meta.k = 20,
                                        clust.seed = 42,
                                        meta.seed = 42)
    dat_train <- flowsom_res[['data']]
    setnames(dat_train, "FlowSOM_metacluster", "FlowSOM_metacluster_original")
    dat_train[, FlowSOM_metacluster := FlowSOM_metacluster_original + meta_multiplication[idx]]
    
    fwrite(dat_train, paste0("Clustered_trainDat_", tp, ".csv"))
    
    # process test data.
    fsom <- flowsom_res[['fsom']]
    
    
    # convert to flowFrame. Pain in the butt
    metadata <- data.frame(name=dimnames(dat_test)[[2]], desc=paste('column',dimnames(dat_test)[[2]],'from dataset'))
    dat_test_ff <- new("flowFrame",
                       exprs=as.matrix(dat_test), # in order to create a flow frame, data needs to be read as matrix
                       parameters=Biobase::AnnotatedDataFrame(metadata))
    
    fsom_2 <- FlowSOM::NewData(fsom, dat_test_ff)
    dat_test <- dat_for_this_tp[dat_for_this_tp$PatientName %in% fold_test$PatientName, ]
    dat_test$FlowSOM_cluster <- fsom_2$map$mapping[,1]
    
    dat_train_code_meta <- data.table(FlowSOM_metacluster=dat_train$FlowSOM_metacluster)
    dat_train_code_meta$FlowSOM_cluster <- fsom$map$mapping[,1]
    dat_train_code_meta <- unique(dat_train_code_meta[, c("FlowSOM_metacluster", "FlowSOM_cluster")])
    dat_test$row_idx <- c(1:nrow(dat_test))
    dat_test <- merge.data.table(dat_test, dat_train_code_meta)
    dat_test <- dat_test[order(row_idx)]
    dat_test$row_idx <- NULL
    fwrite(dat_test, paste0("Clustered_testDat_", tp, ".csv"))
    
  }
  
  if (!file.exists("diffcyt_edgeR_.csv") ||
      !file.exists("diffcyt_limma_.csv")) {
    # Run Diffcyt combining data from all time-points
    dat_train <- lapply(timepoints, function(tp) {
      fread(paste0("Clustered_trainDat_", tp, ".csv"))
    })
    dat_train <- rbindlist(dat_train)
    samp_train <- unique(dat_train$Sample)
    
    markers <- c(lin_markers, func_markers)
    data_ff <- lapply(samp_train, function(samp) {
      input_dat <- dat_train[dat_train$Sample == samp, markers, with = FALSE]
      
      # convert to flowFrame. Pain in the butt
      metadata <- data.frame(name=dimnames(input_dat)[[2]], desc=paste('column',dimnames(input_dat)[[2]],'from dataset'))
      input_dat_ff <- new("flowFrame",
                          exprs=as.matrix(input_dat), # in order to create a flow frame, data needs to be read as matrix
                          parameters=Biobase::AnnotatedDataFrame(metadata))
      return(input_dat_ff)
    })
    names(data_ff) <- samp_train
    
    data_flowset <- as(data_ff, "flowSet")
    rm(data_ff)
    gc()
    
    # setup the experiment info ----
    # we don't want each sample as its own sample as each sample belong to different time-point
    # diffcyt has no awareness of time.
    experiment_info <- unique(dat_train[, c("Sample", "Group", "Timepoint")])
    names(experiment_info) <- c("sample_id", "group_id", "timepoint_id")
    experiment_info[,group_id := ifelse(group_id %in% c("SN", "NR"), 'Non-Reactive', 'Reactive')]
    
    # setup marker info ----
    marker_info <- data.table(marker_name=markers, stringsAsFactors = FALSE)
    marker_info$marker_class <- c(rep("type", length(lin_markers)), rep("state", length(func_markers)))
    
    # Gonna run diffcyt step by step
    d_se <- prepareData(data_flowset, experiment_info, marker_info)
    
    # attach the cluster
    meta_clust <- dat_train$FlowSOM_metacluster
    n_clus <- 80
    meta_clust <- factor(meta_clust, levels = seq_len(n_clus))
    rowData(d_se)$cluster_id <- meta_clust
    
    message("Running edgeR and Limma")
    # Calculate cluster cell counts
    d_counts <- calcCounts(d_se)
    # Calculate cluster medians
    d_medians <- calcMedians(d_se)
    
    # prepare the design matrix, here we're comparing SN against others
    experiment_info$group_id <- factor(experiment_info$group_id, levels = c('Non-Reactive', 'Reactive'))
    design <- createDesignMatrix(experiment_info, cols_design = 'group_id')
    
    contrast <- createContrast(c(0, 1))
    # nrow(contrast) == ncol(design)
    rownames(contrast) <- colnames(design)
    contrast_df <- data.frame(parameters = colnames(design), contrast)
    
    # Run tests
    res_DA <- testDA_edgeR(d_counts, design, contrast, min_cells = 0)
    res_DS <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE, min_cells = 0)
    
    # threshold <- 0.01
    # table(diffcyt::topTable(res_DA, all = TRUE)$p_adj <= threshold)
    # table(diffcyt::topTable(res_DS, all = TRUE)$p_adj <= threshold)
    
    # Export diffcyt result
    top_res_DA <- as.data.frame(diffcyt::topTable(res_DA, all = TRUE))
    fwrite(top_res_DA, "diffcyt_edgeR_.csv")
    
    top_res_DS <- as.data.frame(diffcyt::topTable(res_DS, all = TRUE))
    fwrite(top_res_DS, "diffcyt_limma_.csv")
  }
  
  
}
stopCluster(cl)

# extract features ----
cl <- makeCluster(5) #not to overload your computer
registerDoSNOW(cl)

foreach(i=1:26, .combine = cbind, .packages = c("data.table")) %dopar% {
  source("M:/givanna/cmv/code/extract_diffcyt_features.R")
  timepoints <- c(1:4)
  
  fold_id <- paste0("fold_",i)
  message(paste("Processing", fold_id))
  setwd("M:/givanna/cmv/diffcyt_10x10_20_completeTP_marker_from_LS/")
  setwd(fold_id)
  
  train_dat_list <- list()
  test_dat_list <- list()
  
  # if running diffcyt all tps together
  cnt_dat <- fread("diffcyt_edgeR_.csv")
  marker_dat <- fread("diffcyt_limma_.csv")
  
  cnt_dat <- cnt_dat[order(p_val)]
  marker_dat <- marker_dat[order(p_val)]
  
  cnt_sig <- cnt_dat[1:3,]$cluster_id
  marker_sig <- marker_dat[1:3, c("cluster_id", "marker_id")]
  
  dat_train_list <- lapply(timepoints, function(tp) {
    dat <- fread(paste0("Clustered_trainDat_T", tp, '.csv'))
    dat[, fsom_meta := FlowSOM_metacluster]
    return(dat)
  })
  train_dat <- rbindlist(dat_train_list)
  train_dat_feats <- extract_feats(cnt_sig, marker_sig, train_dat, "fsom_meta")
  
  dat_test_list <- lapply(timepoints, function(tp) {
    dat <- fread(paste0("Clustered_testDat_T", tp, '.csv'))
    dat[, fsom_meta := FlowSOM_metacluster]
    return(dat)
  })
  test_dat <- rbindlist(dat_test_list)
  test_dat_feats <- extract_feats(cnt_sig, marker_sig, test_dat, "fsom_meta")
  
  fwrite(train_dat_feats, "train_data_6imp_feats_diffcyt.csv")
  fwrite(test_dat_feats, "test_data_6imp_feats_diffcyt.csv")
                        
}