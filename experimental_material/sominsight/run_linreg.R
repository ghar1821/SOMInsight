# Select important features from lmms for TrackSOM ----

# This will select top 6, 8, 10, 12 features 
# ready for classification.
# Features must have adjusted p value less than 0.05
library(data.table)
library(stringr)
library(ggpubr)
source("M:/givanna/cmv/code/common_functions.R")

fold_ids <- sapply(c(1:26), function(x) paste0("fold_",x))

tsom_mode <- "10x10_20"
total_n_feat <- 6

main_dir <- paste0("M:/givanna/cmv/tracksom_", tsom_mode, "_completeTP_marker_from_LS")
setwd(main_dir)
save_dir <- paste0("lmms_features_NAmedian_", total_n_feat, "features")
dir.create(save_dir)
# run in serial ----
for (fold_id in fold_ids) {
  # fold_id <- fold_ids[1]
  message(paste("Extracting LMMSDE", total_n_feat, "important features for", fold_id, "of", tsom_mode))
  
  # setwd("M:/givanna/cmv/tracksom_10x10_20_noMerge_completeTP_withCDK1")
  setwd(main_dir)
  setwd(fold_id)
  # read training and test data, and isolate the features ----
  # read training and test data, and isolate the features ----
  dat_train <- lapply(c(1:4), function(i) {
    dat_file <- paste0("Clustered_trainDat_Timepoint_T",i,".csv")
    fread(dat_file)
  })
  dat_train <- rbindlist(dat_train)
  
  dat_test <- lapply(c(1:4), function(i) {
    dat_file <- paste0("Clustered_testDat_Timepoint_T",i,".csv")
    fread(dat_file)
  })
  dat_test <- rbindlist(dat_test)
  
  use_proportion <- TRUE
  
  # setwd("M:/givanna/cmv/tracksom_10x10_20_noMerge_completeTP_withCDK1")
  setwd(paste0(main_dir, "/lmms_features_NAmedian/"))
  setwd(fold_id)
  # TODO use me if you run count and markers independent of one another
  if (use_proportion) 
    abd_dat <- fread("prop_op_all.csv")
  else
    abd_dat <- fread("cnt_op_all.csv")
  
  marker_dat <- fread("marker_op_all.csv")
  
  # 1st, set the key to the adjusted p value column
  setkeyv(abd_dat, c("adj_PVal", "PVal"))
  setkeyv(marker_dat, c("adj_PVal", "PVal"))
  
  # 2nd, remove NA rows
  abd_dat <- na.omit(abd_dat)
  marker_dat <- na.omit(marker_dat)
  
  # sig_feats <- unlist(sig_feats, use.names = FALSE)
  # note adjust me!
  n_feat <- total_n_feat/2
  sig_feats <- c(abd_dat[1:n_feat,]$metacluster, marker_dat[1:n_feat,]$metacluster)
  
  timepoints <- c("T1", "T2", "T3", "T4")
  
  patients_train <- unique(dat_train$PatientName)
  patients_test <- unique(dat_test$PatientName)
  
  
  count_per_sample <- dat_train[, .(count = .N), by=c("Sample")]
  
  corr_scores_train <- extract_dynamic_features(sig_feats = sig_feats,
                                                patients = patients_train,
                                                dat = dat_train,
                                                timepoints = timepoints,
                                                method = 'linear_regression',
                                                use_proportion = use_proportion,
                                                count_per_sample = count_per_sample)
  
  
  count_per_sample <- dat_test[, .(count = .N), by=c("Sample")]
  corr_scores_test <- extract_dynamic_features(sig_feats = sig_feats,
                                               patients = patients_test,
                                               dat = dat_test,
                                               timepoints = timepoints,
                                               method = 'linear_regression',
                                               use_proportion = use_proportion,
                                               count_per_sample = count_per_sample)
  
  # transform to structure understandable by classifier
  corr_scores_train <- dcast(corr_scores_train, PatientName ~ variable)
  corr_scores_test <- dcast(corr_scores_test, PatientName ~ variable)
  
  setwd(main_dir)
  setwd(save_dir)
  dir.create(fold_id)
  setwd(fold_id)
  
  if (use_proportion) {
    fwrite(corr_scores_train, "train_data_imp_feats_prop_linreg.csv")
    fwrite(corr_scores_test, "test_data_imp_feats_prop_linreg.csv")
  } else {
    fwrite(corr_scores_train, "train_data_imp_feats_cnt_linreg.csv")
    fwrite(corr_scores_test, "test_data_imp_feats_cnt_linreg.csv")
  }
}
