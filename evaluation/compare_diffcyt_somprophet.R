library(data.table)
library(Spectre)
library(ggplot2)
library(viridis)
library(scales)

classifiers <- c("DT", "KNN", "LR", "NB", "RF", "SVM")
folds <- c(1:26)

# plot all functional markers ----
setwd(main_dir)
setwd("tracksom_10x10_20_marker_from_LS/lmms_features_NAmedian_6features/classification_res_proportion/")
dat_test_list <- lapply(folds, function(fold) {
  dat_classifiers <- lapply(classifiers, function(cl) {
    d <- fread(paste0("acc_conf_mat_", cl, "_fold_", fold, '.csv'), 
               drop = 1,
               header = TRUE)
    return(d)
  })
  dat_classifiers <- rbindlist(dat_classifiers)
  return(dat_classifiers)
})
dat_test_all_func <- rbindlist(dat_test_list)
dat_test_all_func$marker <- "SOMProphet's temporally dynamic immune features"

setwd(main_dir)
setwd("diffcyt_10x10_20_marker_from_LS/classification_res_allTPTogether/")
dat_test_list <- lapply(folds, function(fold) {
  dat_classifiers <- lapply(classifiers, function(cl) {
    d <- fread(paste0("acc_conf_mat_", cl, "_fold_", fold, '.csv'), 
               drop = 1,
               header = TRUE)
    return(d)
  })
  dat_classifiers <- rbindlist(dat_classifiers)
  return(dat_classifiers)
})
dat_test_diffcyt <- rbindlist(dat_test_list)
dat_test_diffcyt$marker <- "Diffcyt's single static immune features"

dat_test <- rbind(dat_test_all_func, dat_test_diffcyt)
dat_test[, correct := ifelse(Predicted == Reference, "YES", "NO")]
dat_test$correct <- factor(dat_test$correct, levels = c('NO', 'YES'))
dat_test_stat <- dat_test[, .(cnt = .N), by=c('Classifier', 'correct', 'marker')]
dat_test_stat[, proportion := cnt/26]

# replace space with new line for axis
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
dat_test_stat$marker <-
  factor(
    dat_test_stat$marker,
    levels = c(
      "SOMProphet's temporally dynamic immune features",
      "Diffcyt's single static immune features"
    )
  )

setwd(main_dir)
setwd("comparison")

# normal accuracy bar chart.
dat_acc <- dat_test_stat[dat_test_stat$correct == 'YES']
dat_acc[, round_prop := round(proportion, 2)]

ggplot(data=dat_acc, aes(x=marker, y=proportion, fill=marker)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round_prop), vjust=1.6, color="white", size=5) +
  facet_wrap(~ Classifier) +
  scale_x_discrete(labels = addline_format(c("SOMProphet's temporally dynamic immune features",
                                             "Diffcyt's single static immune features"))) +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  labs(y = "Classification accuracy score", x = 'Features') + 
  ggtitle("The classification accuracy of trained classifiers") +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.position = 'none',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12),
        title = element_text(size=14)) +
  scale_y_continuous(breaks=pretty_breaks(n=10), limits = c(0,1))
ggsave("acc_bar.pdf",
       width = 7,
       height = 7)
fwrite(dat_test_stat, "td_vs_pit.csv")