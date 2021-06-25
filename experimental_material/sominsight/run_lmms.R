# devtools::install_github("cran/lmms")

library(lmms)
library(data.table)
library(tools)

source("~/Documents/GitHub/tracksom-differential/v2/common_functions.R")

# func_markers <- c("HCMV_IE", "HCMV_pp65")
func_markers <- c("HCMV_IE", "HCMV_pp65", "CD25", "CD69", "HLA_DR", 'CD38', "CD279_PD1", "Tim3",
                  "CD278", "Granzyme_B", "CD314_NKG2D", "CD86", "CD274")

fold_ids <- sapply(c(1:26), function(x) paste0("fold_",x))
# fold_id <- fold_ids[1]

tsom_mode <- "10x10_20"


for (fold_id in fold_ids) {
  message(paste("Running LMMSDE on", fold_id, "for", tsom_mode))
  setwd(paste0("M:/givanna/cmv/tracksom_", tsom_mode, "_noMerge_completeTP"))
  setwd(fold_id)
  
  message("Reading data")
  # Read clustered train data
  train_dat <- lapply(c(1:4), function(i) {
    dat_file <- paste0("Clustered_trainDat_Timepoint_T",i,".csv")
    return(fread(dat_file))
  })
  train_dat <- rbindlist(train_dat)
  
  # lmms for cluster count ----
  # Unlike other tools, LMMS requires each row to be a sample, and each column to be a gene.
  # Similar to cytometry. So we extract as per normal (row=gene, col=sample), then rotate.
  train_dat_cnt <- extract_cnt_perMetaSample(train_dat, use_proportion = FALSE)
  train_dat_cnt <- t(train_dat_cnt)
  
  train_dat_prop <- extract_cnt_perMetaSample(train_dat, use_proportion = TRUE)
  train_dat_prop <- t(train_dat_prop)
  
  
  train_dat_marker <- extract_markerStats_perMetaSample(dat = train_dat, markers = func_markers, stat_method = 'median')
  train_dat_marker <- t(train_dat_marker)
  
  # Extract timepoint, group for each sample
  samp <- data.table(Sample=rownames(train_dat_marker))
  samp$id <- c(1:nrow(samp))
  samp_dets <- unique(train_dat[, c("Sample", "Group", "Timepoint", "PatientName")])
  samp_dets <- merge(samp, samp_dets)
  samp_dets <- samp_dets[order(id)]
  
  # timepoints have to be numeric vector for lmms
  samp_dets[, Time := as.numeric(gsub("T", "", Timepoint))]
  
  samp_dets[,group_reactive := ifelse(Group %in% c("SN", "NR"), 'Non-Reactive', 'Reactive')]
  samp_dets$Group <- NULL
  
  
  # change the experiment to all (option 1), longitudinal1 (option 3), longitudinal2 (option 4).
  # we'll use cubic for now.
  exp_type <- 'all'
  
  # setup save dir
  save_dir <- paste0("M:/givanna/cmv/tracksom_", tsom_mode, "_noMerge_completeTP/lmms_features_NAmedian/", fold_id)
  
  dir.create(save_dir, recursive = TRUE)
  setwd(save_dir)
  
  message("Running LMMSDE")
  # LMMS for proportion ----
  train_dat_prop_filtered <- remove.low.cv(train_dat_prop)
  
  lmms_out <- lmmsDE(data = train_dat_prop_filtered,
                     time = samp_dets$Time,
                     group = samp_dets$group_reactive,
                     sampleID = samp_dets$Sample,
                     type = 'time*group',
                     experiment = exp_type,
                     basis = 'cubic')
  
  res_filename <- paste0("prop_op_", exp_type,".csv")
  res <- data.table(lmms_out@DE)
  names(res) <- c("metacluster", "PVal", "adj_PVal")
  fwrite(res, res_filename)
  
  # LMMS for count ----
  train_dat_cnt_filtered <- remove.low.cv(train_dat_cnt)
  
  
  lmms_out <- lmmsDE(data = train_dat_cnt_filtered,
                     time = samp_dets$Time,
                     group = samp_dets$group_reactive,
                     sampleID = samp_dets$Sample,
                     type = 'time*group',
                     experiment = exp_type,
                     basis = 'cubic')
  
  res_filename <- paste0("cnt_op_", exp_type,".csv")
  res <- data.table(lmms_out@DE)
  names(res) <- c("metacluster", "PVal", "adj_PVal")
  fwrite(res, res_filename)
  
  # LMMS for markers ----
  train_dat_marker_filtered <- remove.low.cv(train_dat_marker, cutoff = 0.2)
  
  lmms_out <- lmmsDE(data = train_dat_marker_filtered,
                     time = samp_dets$Time,
                     group = samp_dets$group_reactive,
                     sampleID = samp_dets$Sample,
                     type = 'time*group',
                     experiment = exp_type,
                     basis = 'cubic')
  res_filename <- paste0("marker_op_", exp_type,".csv")
  res <- data.table(lmms_out@DE)
  names(res) <- c("metacluster", "PVal", "adj_PVal")
  fwrite(res, res_filename)
}



