library(data.table)

extract_feats <- function(cnt_sig, marker_sig, dat, meta_col_name) {
  
  patients <- unique(dat$PatientName)
  
  # count data
  cnt_dat <- dat[dat[[meta_col_name]] %in% cnt_sig,]
  cnt_dat <- cnt_dat[, .(count = .N), by=c(meta_col_name, "PatientName")]
  cnt_dat <- reshape(cnt_dat,
                           timevar = meta_col_name, 
                           idvar = 'PatientName',
                           direction = 'wide')
  colnames(cnt_dat) <- gsub(".*\\.", paste0(meta_col_name, "_"), colnames(cnt_dat))
  
  # add in missing patients
  cnt_pat <- unique(cnt_dat$PatientName)
  missing_pats <- setdiff(patients, cnt_pat)
  for (p in missing_pats) {
    row <- data.table(PatientName = p)
    for (c in cnt_sig) {
      row[[paste0(meta_col_name, "_", c)]] <- 0
    }
    cnt_dat <- rbind(cnt_dat, row, fill=TRUE)
  }
  
  
  # marker data
  marker_dat <- lapply(c(1:nrow(marker_sig)), function(row_idx) {
    sig_marker <- marker_sig[row_idx,]
    dat <- dat[dat[[meta_col_name]] == sig_marker$cluster_id, 
                     c("PatientName", sig_marker$marker_id), with=FALSE]
    if (nrow(dat) == 0) return(NULL)
    dat <- dat[, lapply(.SD, median), by=c("PatientName"), .SDcols = sig_marker$marker_id]
    setnames(dat, sig_marker$marker_id, paste0(meta_col_name, "_", sig_marker$cluster_id, "_", sig_marker$marker_id))
    return(dat)
  })
  # remove any content that is null. otherwise reduce will return empty data.table
  marker_dat[sapply(marker_dat, is.null)] <- NULL
  
  marker_dat <- Reduce(merge,marker_dat)
  
  # if this above return nothing, it will be null, which cause the merge to fail.
  if (is.null(marker_dat))
    marker_dat <- data.table(PatientName=character())
  
  # add in missing patients
  marker_pat <- unique(marker_dat$PatientName)
  missing_pats <- setdiff(patients, marker_pat)
  for (p in missing_pats) {
    row <- data.table(PatientName = p)
    for (idx in c(1: nrow(marker_sig))) {
      marker_sig_row <- marker_sig[idx]
      row[[paste0(meta_col_name, "_", marker_sig_row$cluster_id, "_", marker_sig_row$marker_id)]] <- NA
    }
    marker_dat <- rbind(marker_dat, row, fill=TRUE)
  }
  
  dat_sig <- merge.data.table(cnt_dat, marker_dat, by='PatientName')
  return(dat_sig)
}





