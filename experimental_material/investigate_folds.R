library(data.table)
library(stringr)
library(ggplot2)

setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/loocv_folds/")

fold_dets <- fread("fold_11.csv")

setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/investigate_marker_from_LS/fold_11/")

dat <- lapply(c(1:4), function(i) fread(paste0("Clustered_trainDat_Timepoint_T", i, ".csv")))

dat <- rbindlist(dat)

setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tracksom_10x10_20_marker_from_LS/lmms_features_NAmedian/fold_11")
abd_dat <- fread("prop_op_all.csv")
marker_dat <- fread("marker_op_all.csv")

# select top 5% of features based on adjusted p value ----
# 1st, set the key to the adjusted p value column
setkeyv(abd_dat, c("adj_PVal", "PVal"))
setkeyv(marker_dat, c("adj_PVal", "PVal"))

# 2nd, remove NA rows
abd_dat <- na.omit(abd_dat)
marker_dat <- na.omit(marker_dat)

# sig_feats <- unlist(sig_feats, use.names = FALSE)
# note adjust me!
total_n_feat <- 6
n_feat <- total_n_feat/2
sig_feats <- c(abd_dat[1:n_feat,]$metacluster, marker_dat[1:n_feat,]$metacluster)

timepoints <- c("T1", "T2", "T3", "T4")

count_per_sample <- dat[, .(count = .N), by=c("Sample")]
sample_patient_timepoint_mapping <- unique(dat[, c("Sample", "PatientName", "Timepoint")])

patients <- unique(dat$PatientName)
use_proportion <- TRUE

dat_sig_feat_list <- lapply(sig_feats, function(sig_feat) {
  # this is cluster count
  if (isFALSE(grepl("_", sig_feat))) {
    corr_scores <- lapply(patients, function(pat) {
      pat_dat <- dat[dat$PatientName == pat &
                       TrackSOM_metacluster_lineage_tracking == sig_feat]
      count_per_tp <- pat_dat[, .(count = .N), by = c("Timepoint")]
      
      # fill in missing time-points
      missing_tp <- setdiff(timepoints, count_per_tp$Timepoint)
      missing_dat <- data.table(Timepoint=missing_tp, count=rep(0, length(missing_tp)))
      count_per_tp <- rbind(count_per_tp, missing_dat)
      
      if (use_proportion) {
        sample_timepoint_mapping <- sample_patient_timepoint_mapping[sample_patient_timepoint_mapping$PatientName == pat]
        
        cell_proportions <- sapply(c(1: nrow(count_per_tp)), function(row_idx) {
          row <- count_per_tp[row_idx, ]
          cnt <- row$count
          
          sample_id <- sample_timepoint_mapping[sample_timepoint_mapping$Timepoint == row$Timepoint]$Sample
          sample_size <- count_per_sample[count_per_sample$Sample == sample_id,]$count
          return(cnt/sample_size)
        })
        count_per_tp$count <- cell_proportions
      }
      
      # substitute timepoint with numeric equivalent
      count_per_tp$Timepoint <- match(count_per_tp$Timepoint, timepoints)
      count_per_tp$PatientName <- pat
      count_per_tp$feat <- sig_feat
      
      setnames(count_per_tp, 'count', 'value')
      return(count_per_tp)
    })
  } else {
    feat_split <- str_split(sig_feat, "_")[[1]]
    last_idx <- length(feat_split) - 1
    marker <- paste(feat_split[c(1: last_idx)], collapse = "_")
    cluster <- feat_split[length(feat_split)]
    
    corr_scores <- lapply(patients, function(pat) {
      pat_dat <- dat[dat$PatientName == pat &
                       TrackSOM_metacluster_lineage_tracking == cluster]
      median_per_tp <- pat_dat[, lapply(.SD, median), .SDcols=c(marker), by="Timepoint"]
      
      # fill in missing time-points
      missing_tp <- setdiff(timepoints, median_per_tp$Timepoint)
      missing_dat <- data.table(Timepoint=missing_tp)
      missing_dat[[marker]] <- rep(0, length(missing_tp))
      median_per_tp <- rbind(median_per_tp, missing_dat)
      
      # substitute timepoint with numeric equivalent
      median_per_tp$Timepoint <- match(median_per_tp$Timepoint, timepoints)
      median_per_tp$PatientName <- pat
      median_per_tp$feat <- sig_feat
      setnames(median_per_tp, marker, 'value')
      
      return(median_per_tp)
    })
  }
  corr_scores <- rbindlist(corr_scores)
  return(corr_scores)
})
dat_sig_feat <- rbindlist(dat_sig_feat_list)

pat_grp <- unique(dat[, c("Group", "PatientName")])

unique_feats <- unique(dat_sig_feat$feat)
dat_sig_marker_feats <- dat_sig_feat[dat_sig_feat$feat %in% unique_feats[4:6]]
dat_sig_marker_feats <- merge(dat_sig_marker_feats, pat_grp, by='PatientName')

setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/investigate_marker_from_LS/fold_11/")

fwrite(dat_sig_feat, "significant_immune_features.csv")

ggplot(data=dat_sig_marker_feats, aes(x=Timepoint, y=value, color=Group)) +
  geom_line() +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) + 
  facet_wrap(feat ~ PatientName)
ggsave("imp_markers.png", p, height=20, width=20)

dat_sig_cnt_feats <- dat_sig_feat[dat_sig_feat$feat %in% unique_feats[1:3]]
dat_sig_cnt_feats <- merge(dat_sig_cnt_feats, pat_grp, by='PatientName')

ggplot(data=dat_sig_cnt_feats, aes(x=Timepoint, y=value, color=Group)) +
  geom_line() +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) + 
  facet_wrap(feat ~ PatientName)
ggsave("imp_cnt.png", p, height=20, width=20)

# compute residuals
corr_scores <- lapply(sig_feats, function(sig_feat) {
  print(sig_feat)
  # this is cluster count
  if (isFALSE(grepl("_", sig_feat))) {
    corr_scores <- lapply(patients, function(pat) {
      pat_dat <- dat[dat$PatientName == pat &
                       TrackSOM_metacluster_lineage_tracking == sig_feat]
      count_per_tp <- pat_dat[, .(count = .N), by = c("Timepoint")]
      
      # fill in missing time-points
      missing_tp <- setdiff(timepoints, count_per_tp$Timepoint)
      missing_dat <- data.table(Timepoint=missing_tp, count=rep(0, length(missing_tp)))
      count_per_tp <- rbind(count_per_tp, missing_dat)
      
      if (use_proportion) {
        sample_timepoint_mapping <- sample_patient_timepoint_mapping[sample_patient_timepoint_mapping$PatientName == pat]
        
        cell_proportions <- sapply(c(1: nrow(count_per_tp)), function(row_idx) {
          row <- count_per_tp[row_idx, ]
          cnt <- row$count
          
          sample_id <- sample_timepoint_mapping[sample_timepoint_mapping$Timepoint == row$Timepoint]$Sample
          sample_size <- count_per_sample[count_per_sample$Sample == sample_id,]$count
          return(cnt/sample_size)
        })
        count_per_tp$count <- cell_proportions
      }
      
      # substitute timepoint with numeric equivalent
      count_per_tp$Timepoint <- match(count_per_tp$Timepoint, timepoints)
      corr_score <- lm(formula = "count ~ Timepoint", data = count_per_tp)
      corr_score <- resid(corr_score)
      corr_score <- data.table(residual=corr_score, timepoints=c(1:4))
      corr_score$PatientName <- pat
      corr_score$feat <- sig_feat
      return(corr_score)
    })
  } else {
    feat_split <- str_split(sig_feat, "_")[[1]]
    last_idx <- length(feat_split) - 1
    marker <- paste(feat_split[c(1: last_idx)], collapse = "_")
    cluster <- feat_split[length(feat_split)]
    
    corr_scores <- lapply(patients, function(pat) {
      pat_dat <- dat[dat$PatientName == pat &
                       TrackSOM_metacluster_lineage_tracking == cluster]
      median_per_tp <- pat_dat[, lapply(.SD, median), .SDcols=c(marker), by="Timepoint"]
      
      # fill in missing time-points
      missing_tp <- setdiff(timepoints, median_per_tp$Timepoint)
      missing_dat <- data.table(Timepoint=missing_tp)
      missing_dat[[marker]] <- rep(0, length(missing_tp))
      median_per_tp <- rbind(median_per_tp, missing_dat)
      
      # substitute timepoint with numeric equivalent
      median_per_tp$Timepoint <- match(median_per_tp$Timepoint, timepoints)
      
      corr_score <- lm(formula = paste0(marker, " ~ Timepoint"), data = median_per_tp)
      corr_score <- resid(corr_score)
      corr_score <- data.table(residual=corr_score, timepoints=c(1:4))
      corr_score$PatientName <- pat
      corr_score$feat <- sig_feat
      
      return(corr_score)
    })
  }
  
  corr_scores_dt <-rbindlist(corr_scores)
  return(corr_scores_dt)
})

corr_scores_list <- corr_scores
corr_scores <- rbindlist(corr_scores_list)

feats <- unique(corr_scores$feat)

titles <- sapply(feats, function(f) {
  if (grepl("_", f)) {
    marker_and_clust <- stringr::str_split(f, "_")[[1]]
    num_elem <- length(marker_and_clust)
    clust <- marker_and_clust[num_elem]
    marker <- paste(marker_and_clust[1: (num_elem-1)], collapse = "_")
    return(paste("Median expression of", marker, "in meta-cluster", clust))
  } else {
    return(paste("Proportion of cells in meta-cluster", f))
  }
})

fwrite(corr_scores, "residual.csv")

for (f in feats) {
  corr_score <- corr_scores[corr_scores$feat == f,]
  
  p <- ggplot(corr_score, aes(x=timepoints, y=residual)) +
    geom_point() +
    geom_line() + 
    facet_wrap(~ PatientName) +
    ggtitle(titles[[f]]) +
    labs(y = 'Residuals', x = 'Time-points') +
    theme_clean()
  ggsave(paste0("residual_", f ,".pdf"), p)
}


