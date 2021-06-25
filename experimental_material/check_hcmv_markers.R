library(data.table)
library(stringr)
library(ggplot2)

# check the positions of the HCMV markers ----
dat_dir <- "~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tracksom_10x10_20_marker_from_LS/lmms_features_NAmedian/"
features <- lapply(folds, function(fold) {
  setwd(dat_dir)
  setwd(paste0("fold_", fold))
  dat <- fread("marker_op_all.csv")
  dat <- dat[order(PVal)]
  dat$rank <- c(1:nrow(dat))
  # print(nrow(dat))
  dat_hcmv <- dat[grepl("HCMV", dat$metacluster, fixed=TRUE),]
  dat_hcmv$fold <- fold
  return(dat_hcmv)
})
features <- rbindlist(features)
features[, .(cnt = .N), by='fold']

feat_pp65 <- features[grepl("pp65", features$metacluster, fixed=TRUE),]
feat_pp65$marker <- "HCMV pp65"
mean(feat_pp65$rank)
min(feat_pp65$rank)
feat_ie <- features[grepl("IE", features$metacluster, fixed=TRUE),]
feat_ie$marker <- "HCMV IE"
mean(feat_ie$rank)
min(feat_ie$rank)

# check ranking within top 10
nrow(feat_pp65[feat_pp65$rank <= 10,])/nrow(feat_pp65)
nrow(feat_ie[feat_ie$rank <= 10,])/nrow(feat_ie)

# plot cdf for both
feat_hcmv <- rbind(feat_pp65, feat_ie)
ggplot(feat_hcmv, aes(x = rank, colour = marker)) + 
  stat_ecdf(geom = "step") +
  theme_bw() +
  labs(y = "Cumulative proportion", x = 'Ranking', colour = 'Marker') +
  ggtitle("Cumulative proportion of HCMV-specific temporal state features ranking") +
  scale_color_brewer(palette="Set1") +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.position = 'bottom',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12),
        title = element_text(size=12)) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_continuous(breaks=pretty_breaks(n=10)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
setwd("~/Dropbox (Sydney Uni)/tracksom_differential/")
ggsave("hcmv_marker_ranking.pdf",
       width = 7,
       height = 7)

ggplot(feat_hcmv, aes(x = PVal, colour = marker)) + 
  stat_ecdf(geom = "step") +
  theme_bw() +
  labs(y = "Cumulative proportion", x = 'P-value', colour = 'Marker') +
  ggtitle("Cumulative proportion of HCMV-specific temporal state features' P-value") +
  scale_color_brewer(palette="Set1") +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.position = 'bottom',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12),
        title = element_text(size=12)) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_continuous(breaks=pretty_breaks(n=10)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
setwd("~/Dropbox (Sydney Uni)/tracksom_differential/")
ggsave("hcmv_marker_pval.pdf",
       width = 7,
       height = 7)
