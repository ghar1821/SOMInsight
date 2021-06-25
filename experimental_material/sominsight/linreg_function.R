#' Run linear regression
#' 
#' @param top_n_feats number of features to extract, ordered based on p-value. Used if pval_thres is NULL.
#'  Applicable to both proportion and state data.
#' @param pval_thres p-value threshold which determine which feature to extract. Used if top_n_feats is NULL.
run_linreg <- function(dat, timepoints,
                       top_n_feats = NULL, pval_thres = NULL,
                       use_proportion = TRUE) {
  if (use_proportion) {
    abd_dat <- fread("prop_op_all.csv")
  }
  else {
    abd_dat <- fread("cnt_op_all.csv")
  }
    
  marker_dat <- fread("marker_op_all.csv")
  
  setkeyv(abd_dat, c("adj_PVal", "PVal"))
  setkeyv(marker_dat, c("adj_PVal", "PVal"))
  
  # 2nd, remove NA rows
  abd_dat <- na.omit(abd_dat)
  marker_dat <- na.omit(marker_dat)
  
  if (is.null(top_n_feats) & !is.null(pval_thres)) {
    # extract based on p_value
    abd_feats <- abd_dat[abd_dat$PVal < pval_thres]
    marker_feats <- marker_dat[marker_dat$PVal < pval_thres]
    sig_feats <- c(abd_feats$metacluster, marker_feats$metacluster)
  } else if (!is.null(top_n_feats) & is.null(pval_thres)) {
    sig_feats <- c(abd_dat[1:top_n_feats,]$metacluster, marker_dat[1:top_n_feats,]$metacluster)
  }
  
  patients <- unique(dat$PatientName)
  
  count_per_sample <- dat[, .(count = .N), by=c("Sample")]
  
  corr_scores <- extract_dynamic_features(sig_feats = sig_feats,
                                          patients = patients,
                                          dat = dat,
                                          timepoints = timepoints,
                                          method = 'linear_regression',
                                          use_proportion = use_proportion,
                                          count_per_sample = count_per_sample)
  
  
  # transform to structure understandable by classifier
  corr_scores <- dcast(corr_scores, PatientName ~ variable)
  
  if (use_proportion) {
    fwrite(corr_scores, "dat_prop_linreg.csv")
  } else {
    fwrite(corr_scores, "dat_cnt_linreg.csv")
  }
}
