# All the common operations required among all approaches

extract_cnt_perMetaSample <- function(dat, use_proportion=TRUE) {
  meta_count_per_sample <- dat[, .(count = .N), by=c("Sample", "TrackSOM_metacluster_lineage_tracking")]
  
  if (use_proportion) {
    # compute proportion over samples rather than over patient or timepoint.
    # if over time-point, patient may have low proportion just because he donated less cells.
    # if over patients, patient may have low proportion just because that time-point there are much less cells.
    count_per_samp <- dat[, .(count = .N), by=c("Sample")]
    
    meta_prop_per_sample <- sapply(c(1: nrow(meta_count_per_sample)), function(row_idx) {
      row <- meta_count_per_sample[row_idx,]
      cnt_for_samp <- count_per_samp[count_per_samp$Sample == row$Sample]$count
      proportion <- row$count / cnt_for_samp
      return(proportion)
    })
    meta_count_per_sample$count <- meta_prop_per_sample
  }
  
  meta_count_mat <- reshape(meta_count_per_sample, 
                            timevar = "Sample", 
                            idvar = 'TrackSOM_metacluster_lineage_tracking',
                            direction = 'wide')
  colnames(meta_count_mat) = gsub(".*\\.", "", colnames(meta_count_mat))
  
  # store for now
  meta_clust <- meta_count_mat$TrackSOM_metacluster_lineage_tracking
  
  meta_count_mat[, TrackSOM_metacluster_lineage_tracking := NULL]
  
  # replace NA with 0
  meta_count_mat[is.na(meta_count_mat)] <- 0
  
  meta_count_mat <- as.matrix(meta_count_mat)
  rownames(meta_count_mat) <- meta_clust
  
  return(meta_count_mat)
}

extract_markerStats_perMetaSample <- function(dat, markers, stat_method = c("median", "mean", "sd")) {
  
  stat_method <- match.arg(stat_method)
  
  if (stat_method == 'median') {
    meta_med_marker_samp <- dat[, lapply(.SD, median), by=c("Sample", "TrackSOM_metacluster_lineage_tracking"), .SDcols = markers]
  } else if (stat_method == 'mean') {
    meta_med_marker_samp <- dat[, lapply(.SD, mean), by=c("Sample", "TrackSOM_metacluster_lineage_tracking"), .SDcols = markers]
  } else if (stat_method == 'sd') {
    meta_med_marker_samp <- dat[, lapply(.SD, sd), by=c("Sample", "TrackSOM_metacluster_lineage_tracking"), .SDcols = markers]
  } else {
    stop("stat_method not understood. Must be median or mean.")
  }
  
  # we need to create a data frame where each row is median marker exp of a meta cluster
  meta_med_markers_mats <- lapply(markers, function(marker) {
    cols <- c('Sample', 'TrackSOM_metacluster_lineage_tracking', marker)
    meta_med_one_marker <- meta_med_marker_samp[, cols, with=FALSE]
    
    meta_med_marker_mat <- reshape(meta_med_one_marker, 
                                   timevar = "Sample", 
                                   idvar = 'TrackSOM_metacluster_lineage_tracking',
                                   direction = 'wide')
    colnames(meta_med_marker_mat) = gsub(".*\\.", "", colnames(meta_med_marker_mat))
    
    # concatenate meta cluster and marker
    meta_med_marker_mat[, marker_meta:=paste(marker, TrackSOM_metacluster_lineage_tracking, sep="_")]
    meta_med_marker_mat[, TrackSOM_metacluster_lineage_tracking:=NULL]
    return(meta_med_marker_mat)
  })
  meta_med_markers_mats <- rbindlist(meta_med_markers_mats)
  
  # replace NA with 0
  # meta_med_markers_mats[is.na(meta_med_markers_mats)] <- 0
  row_names <- meta_med_markers_mats$marker_meta
  meta_med_markers_mats[, marker_meta := NULL]
  meta_med_markers_mtx <- as.matrix(meta_med_markers_mats)
  rownames(meta_med_markers_mtx) <- row_names
  
  return(meta_med_markers_mtx)
}

#' Remove low CV
#' 
#' In a longitudinal context, one can be interested only in features that vary over time and 
#' filter out molecules with a low variation coefficient.
#' To do so, we can first naively set a threshold on the variation coefficient and 
#' keep those features that exceed the threshold.
#' 
#' @param X data.table Data to be filtered
#' @param cutoff numeric Default 0.5, the CV cutoff value 
#' 
#' @return filtered data frame
remove.low.cv <- function(X, cutoff = 0.5){
  # X_scaled <- scales::rescale(X, to=c(0,1))
  
  # var.coef
  cv <- unlist(lapply(as.data.frame(X), function(y) {
    return(abs(sd(y, na.rm = TRUE)/mean(y, na.rm=TRUE)))
  }))
  
  
  return(X[,!is.na(cv) & cv > cutoff])
}

remove.low.var <- function(X, cutoff = 0.5){
  # var.coef
  variances <- unlist(lapply(as.data.frame(X), 
                             function(x) abs(var(x))))
  
  return(X[,variances > cutoff])
}

extract_dynamic_features <- function(sig_feats, patients, dat, timepoints,
                                     method = c("correlation", "linear_regression"),
                                     use_proportion = FALSE, count_per_sample = NULL) {
  method <- match.arg(method)
  
  if (use_proportion) {
    sample_patient_timepoint_mapping <- unique(dat[, c("Sample", "PatientName", "Timepoint")])
  }
  
  message(paste("Extracting", paste0(sig_feats, collapse = ",")))
  corr_scores <- lapply(sig_feats, function(sig_feat) {
    # this is cluster count
    message(paste("Processing", sig_feat))
    if (isFALSE(grepl("_", sig_feat))) {
      corr_scores <- sapply(patients, function(pat) {
        message(paste("Processing", pat))
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
        
        if (method == 'linear_regression') {
          corr_score <- lm(formula = "count ~ Timepoint", data = count_per_tp)
          corr_score <- corr_score$coefficients[2]    
        } else {
          corr_score <- cor(count_per_tp$Timepoint, count_per_tp$count, method='spearman')
        }
        return(corr_score)
      })
    } else {
      feat_split <- str_split(sig_feat, "_")[[1]]
      last_idx <- length(feat_split) - 1
      marker <- paste(feat_split[c(1: last_idx)], collapse = "_")
      cluster <- feat_split[length(feat_split)]
      
      corr_scores <- sapply(patients, function(pat) {
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
        
        if (method == 'linear_regression') {
          corr_score <- lm(formula = paste0(marker, " ~ Timepoint"), data = median_per_tp)
          corr_score <- corr_score$coefficients[2]    
        } else {
          corr_score <- cor(median_per_tp$Timepoint, median_per_tp[[marker]], method='spearman')
        }
        return(corr_score)
      })
    }
    
    corr_scores_dt <- data.table(PatientName=patients)
    corr_scores_dt[[sig_feat]] <- corr_scores
    corr_scores_dt <- melt(corr_scores_dt, id.vars="PatientName")
    return(corr_scores_dt)
  })
  corr_scores <- rbindlist(corr_scores)
  return(corr_scores)
}

