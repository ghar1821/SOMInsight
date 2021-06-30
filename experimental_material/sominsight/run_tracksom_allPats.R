# TrackSOM is ran on training data.
# Testing data is mapped onto FlowSOM's SOM nodes.
# Based on the node membership, determine the cell's meta clusters
# Import libraries
library(TrackSOM)
library(Spectre)
library(tools)

# run TrackSOM ----

source('M:/givanna/cmv/code/TrackSOM/TrackSOM.R')
source('M:/givanna/cmv/code/TrackSOM/TrackingFunctions.R')
source('M:/givanna/cmv/code/TrackSOM/UtilityFunctions.R')
source('M:/givanna/cmv/code/TrackSOM/TrackingNoMerging.R')

data_dir <- "M:/givanna/cmv/data_complete_timepoints_only/"
data.files <- list.files(data_dir, ".csv")
dat_per_tp <- lapply(data.files, function(f) fread(paste0(data_dir, f)))
markers <- c("CDK1","CD57", "CD274", "CD159c_NKG2C", "CD161", "FOXP3", "CD34", "CD38", "HCMV_pp65", "CD3", "CD278", "CD304", "CD163", "CD314_NKG2D", "CD86", "CD33", "CD14", "CD127", "CD11C", "CD279_PD1", "CD56", "CD19", "CD45RA", "CD69", "CD4", "CD8A", "CD16", "TCRgd", "CD27", "CD20", "Tim3", "CD45RO", "HCMV_IE", "CD25", "Granzyme_B", "CD197_CCR7", "TCRva7.2", "HLA_DR")
  
message("Running TrackSOM")
tracksom.result <- TrackSOM(inputFiles = dat_per_tp,
                            colsToUse = markers,
                            tracking = TRUE,
                            noMerge = TRUE,  # TODO change me if you want to allow merging
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            nClus = 20,
                            scale = TRUE,
                            dataFileType = "data.frame"  # TODO change me according to file type you have
)
dat_train <- rbindlist(dat_per_tp)
dat_train <- ConcatenateClusteringDetailss(tracksom.result = tracksom.result,
                                         dat = dat_train,
                                         timepoint.col = 'Timepoint',
                                         timepoints = c("T1", "T2", "T3", "T4"))
message("Saving training data")
save_dir <- "M:/givanna/cmv/tsom_all_pats/"
dir.create(save_dir, recursive = TRUE)
setwd(save_dir)
Spectre::write.files(dat_train, file.prefix = 'Clustered_trainDat', divide.by = 'Timepoint')

# draw up some plots ----
library(Spectre)
library(data.table)

# source("M:/givanna/cmv/code/TrackSOM/viz/get.transitions.R")
# source("M:/givanna/cmv/code/TrackSOM/viz/make.network.plot.R")
source("~/Documents/GitHub/FlowSOM-tracking/TrackSOM/viz/get.transitions.R")
source("~/Documents/GitHub/FlowSOM-tracking/TrackSOM/viz/make.network.plot.R")

timepoints <- c("T1", "T2", "T3", "T4")

dat_dir <- "~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tsom_all_pats/"
setwd(dat_dir)

dat <- lapply(timepoints, function(t) {
  return(fread(paste0("Clustered__Timepoint_", t, ".csv")))
})
dat <- rbindlist(dat)

# Remove columns we don't need to save RAM
cols_to_keep <- names(dat)[c(1:38, 42, 43, 44, 45, 52)]
dat <- dat[, cols_to_keep, with=FALSE]
markers <- names(dat)[1:38]

dir.create("network_plots/")
setwd("network_plots/")
make.network.plot(dat = dat,
                  timepoint.col = 'Timepoint',
                  timepoints = timepoints,
                  cluster.col = "TrackSOM_metacluster_lineage_tracking",
                  marker.cols = markers,
                  calculation.type = 'median',
                  graph.layout = 'gem',
                  no_merge = TRUE,
                  file.format = 'png',
                  load.temp.data = FALSE,
                  node.size = 5,
                  arrow.head.gap = 2,
                  arrow.length = 2,
                  standard.colours = 'spectral',
                  min.node.size = 1,
                  max.node.size = 10)


# Infer temporal-dynamic features ----
# source("M:/givanna/cmv/code/common_functions.R")
source("~/Documents/GitHub/tracksom-differential/v2/common_functions.R")
func_markers <- c("HCMV_IE", "HCMV_pp65", "CD25", "CD69", "HLA_DR", 'CD38', "CD279_PD1", "Tim3",
                  "CD278", "Granzyme_B", "CD314_NKG2D", "CD86", "CD274")

train_dat_prop <- extract_cnt_perMetaSample(dat, use_proportion = TRUE)
train_dat_prop <- t(train_dat_prop)

# lmms for marker intensity ----
# markers <- names(dat)[c(1:38)]

train_dat_marker <- extract_markerStats_perMetaSample(dat = dat, markers = func_markers, stat_method = 'median')
train_dat_marker <- t(train_dat_marker)

# Extract timepoint, group for each sample
samp <- data.table(Sample=rownames(train_dat_prop))
samp$id <- c(1:nrow(samp))
samp_dets <- unique(dat_train[, c("Sample", "Group", "Timepoint", "PatientName")])
samp_dets <- merge(samp, samp_dets)
samp_dets <- samp_dets[order(id)]

# timepoints have to be numeric vector for lmms
samp_dets[, Time := as.numeric(gsub("T", "", Timepoint))]

samp_dets[,group_reactive := ifelse(Group %in% c("SN", "NR"), 'Non-Reactivated', 'Reactivated')]
samp_dets$Group <- NULL


# change the experiment to all (option 1), longitudinal1 (option 3), longitudinal2 (option 4).
# we'll use cubic for now.
exp_type <- 'all'
save_dir <- "lmms"
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

# run linear regression for some features ----
setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tsom_all_pats/")
dat <- lapply(c(1:4), function(t) {
  fread(paste0("Clustered__Timepoint_T", t, ".csv"))
})
dat <- rbindlist(dat)
timepoints <- c("T1", "T2", "T3", "T4")
setwd("lmms/")
run_linreg(dat = dat,
           timepoints = timepoints,
           pval_thres = 0.05)

# compare the slope of linear regressions ----
dat_group <- fread("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/patient_groups.csv")
setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tsom_all_pats/lmms/")
dat_linreg <- fread("dat_prop_linreg.csv")

feats <- c("T", "H", "M", "D", "E", "Tim3_T", "Granzyme_B_L")

dat_linreg <- merge(dat_linreg, dat_group, by='PatientName')
dat_linreg <- melt(dat_linreg, id.vars = c("PatientName", "group_reactive"), 
     measure.vars = feats)
dat_linreg_prop <- dat_linreg[dat_linreg$variable %in% c("T", "H", "M", "D", "E")]
dat_linreg_marker <- dat_linreg[dat_linreg$variable %in% c("Tim3_T", "Granzyme_B_L")]

# ks test if required
ks_prop <- lapply(feats, function(feat) {
  
  if (feat %in% feats[6:7]) {
    dat_to_process <- dat_linreg_marker
  } else {
    dat_to_process <- dat_linreg_prop
  }
  
  reactivated <- dat_to_process[dat_to_process$group_reactive == 'Reactivated' &
                                  dat_to_process$variable == feat,]
  nonreactivated <- dat_to_process[dat_to_process$group_reactive == 'Non-Reactivated' &
                                     dat_to_process$variable == feat]
  ks_val <- ks.test(reactivated$value, nonreactivated$value)
  ks_val_dt <- data.table(feat=feat,
                          dval=ks_val$statistic,
                          pval=ks_val$p.value,
                          pval_sig=ks_val$p.value<0.05)
  return(ks_val_dt)
})
ks_prop <- rbindlist(ks_prop)

fwrite(ks_prop, "ks_test.csv")

# plots ----
library(ggplot2)
mapping <- fread("meta_mapping.csv")
dat_linreg_prop <- merge(dat_linreg_prop, mapping, by.x = 'variable',
                         by.y = 'metacluster')

dat_linreg_prop$group_reactive <- factor(dat_linreg_prop$group_reactive, 
                                         levels = c("Reactivated", "Non-Reactivated"))

# cumulative distributions ----
ggplot(dat_linreg_prop,
       aes_string(x = "value", colour = 'group_reactive')) +
  stat_ecdf(size = 1) +
  facet_wrap( ~ population, scales = "free") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_bw() +
  labs(x = 'Temporal proportion feature value',
       y = 'Cumulative proportion',
       colour = 'Class') +
  ggtitle("Cumulative proportion of the temporal proportion features values") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_brewer(palette="Set1") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    title = element_text(size = 16),
    legend.position = 'bottom',
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

ggsave("cum_temp_prop_pop.pdf",
       width = 12,
       height = 8,
       limitsize = FALSE)

dat_linreg_marker$group_reactive <- factor(dat_linreg_marker$group_reactive, 
                                         levels = c("Reactivated", "Non-Reactivated"))
dat_linreg_marker <- merge(dat_linreg_marker, mapping, by.x = 'variable',
                         by.y = 'metacluster')
ggplot(dat_linreg_marker, aes_string(x="value", colour='group_reactive')) +
  stat_ecdf(size=1) +
  facet_wrap(~ population, scales = "free") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_bw() +
  labs(x = 'Temporal state feature value',
       y = 'Cumulative proportion',
       colour = 'Group') +
  ggtitle("Cumulative proportion of the temporal state features values") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_brewer(palette="Set1") +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        title = element_text(size=16),
        legend.position = 'bottom',
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        strip.text = element_text(size=16))

ggsave("cum_temp_state.pdf",
       width = 8,
       height = 5,
       limitsize = FALSE)

# autograph ----
dat_dir <- "~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tsom_all_pats/"
setwd(dat_dir)
timepoints <- c("T1", "T2", "T3", "T4")

dat <- lapply(timepoints, function(t) {
  return(fread(paste0("Clustered__Timepoint_", t, ".csv")))
})
dat <- rbindlist(dat)

# read sample details
# Extract timepoint, group for each sample
samp_dets <- unique(dat[, c("Sample", "Group", "Timepoint", "PatientName")])
samp_dets[,group_reactive := ifelse(Group %in% c("SN", "NR"), 'Non-Reactivated', 'Reactivated')]
samp_dets$Group <- NULL

source("~/Documents/GitHub/tracksom-differential/v2/common_functions.R")
func_markers <- c("HCMV_IE", "HCMV_pp65", "CD25", "CD69", "HLA_DR", 'CD38', "CD279_PD1", "Tim3",
                  "CD278", "Granzyme_B", "CD314_NKG2D", "CD86", "CD274")
train_dat_prop <- extract_cnt_perMetaSample(dat, use_proportion = TRUE)

train_dat_prop_dt <- data.table(train_dat_prop)
train_dat_prop_dt$meta_cluster <- rownames(train_dat_prop)

setwd("lmms")
mapping <- fread("meta_mapping.csv")
meta_cnt <- c("T", "H", "D")
train_dat_prop_dt <- train_dat_prop_dt[train_dat_prop_dt$meta_cluster %in% meta_cnt]
train_dat_prop_dt <-
  melt(
    train_dat_prop_dt,
    id.vars = "meta_cluster",
    variable.name = 'Sample',
    value.name = 'proportion'
  )
train_dat_prop_dt <- merge(train_dat_prop_dt, samp_dets)

setwd("..")
library(RColorBrewer)
# NK cells
train_dat_prop_dt_pop <- train_dat_prop_dt[train_dat_prop_dt$meta_cluster == 'T']
p <- make.autograph(dat = train_dat_prop_dt_pop,
               x.axis = "Timepoint",
               y.axis = 'proportion',
               colours = c('blue', 'red'),
               colour.by = 'group_reactive',
               violin = FALSE,
               title = "Proportion of NK cells",
               y.axis.label = 'Proportion of cells',
               x.axis.label = 'Time-point',
               filename = "nk_cells.pdf",
               dot.size = 4
) + facet_grid(. ~ group_reactive) +
  theme(strip.text = element_text(colour="black", size=12, face = 'bold'),
        legend.position = 'none')
ggsave(plot = p, filename = "nk_cells.pdf",
       width = 7, height = 5)

# CD8 T cells
train_dat_prop_dt_pop <- train_dat_prop_dt[train_dat_prop_dt$meta_cluster == 'H']
p <- make.autograph(dat = train_dat_prop_dt_pop,
                    x.axis = "Timepoint",
                    y.axis = 'proportion',
                    colours = c('blue', 'red'),
                    colour.by = 'group_reactive',
                    violin = FALSE,
                    title = "Proportion of CD8+ T cells",
                    y.axis.label = 'Proportion of cells',
                    x.axis.label = 'Time-point',
                    filename = "cd8_tcells.pdf",
                    dot.size = 4
) + facet_grid(. ~ group_reactive) +
  theme(strip.text = element_text(colour="black", size=12, face = 'bold'),
        legend.position = 'none')
ggsave(plot = p, filename = "cd8_tcells.pdf",
       width = 7, height = 5)

# Monocytes
train_dat_prop_dt_pop <- train_dat_prop_dt[train_dat_prop_dt$meta_cluster == 'D']
p <- make.autograph(dat = train_dat_prop_dt_pop,
                    x.axis = "Timepoint",
                    y.axis = 'proportion',
                    colours = c('blue', 'red'),
                    colour.by = 'group_reactive',
                    violin = FALSE,
                    title = "Proportion of Classical Monocytes",
                    y.axis.label = 'Proportion of cells',
                    x.axis.label = 'Time-point',
                    filename = "monocytes.pdf",
                    dot.size = 4
) + facet_grid(. ~ group_reactive) +
  theme(strip.text = element_text(colour="black", size=12, face = 'bold'),
        legend.position = 'none')
ggsave(plot = p, filename = "monocytes.pdf",
       width = 7, height = 5)

# Markers 
train_dat_marker <- extract_markerStats_perMetaSample(dat = dat, markers = func_markers, stat_method = 'median')

train_dat_marker_dt <- data.table(train_dat_marker)
train_dat_marker_dt$meta_cluster <- rownames(train_dat_marker)

mapping <- fread("meta_mapping.csv")
meta_marker <- c("Tim3_T")
train_dat_marker_dt <- train_dat_marker_dt[train_dat_marker_dt$meta_cluster %in% meta_marker]
train_dat_marker_dt <-
  melt(
    train_dat_marker_dt,
    id.vars = "meta_cluster",
    variable.name = 'Sample',
    value.name = 'median_marker'
  )
train_dat_marker_dt <- merge(train_dat_marker_dt, samp_dets)
p <- make.autograph(dat = train_dat_marker_dt,
                    x.axis = "Timepoint",
                    y.axis = 'median_marker',
                    colours = c('blue', 'red'),
                    colour.by = 'group_reactive',
                    violin = FALSE,
                    title = "Median expression of TIM-3 marker in NK cells",
                    y.axis.label = 'Median expression of TIM-3 marker',
                    x.axis.label = 'Time-point',
                    filename = "tim3_nk.pdf",
                    dot.size = 4
) + facet_grid(. ~ group_reactive) +
  theme(strip.text = element_text(colour="black", size=12, face = 'bold'),
        legend.position = 'none')
ggsave(plot = p, filename = "tim3_nk.pdf",
       width = 7, height = 5)
  

# investigate meta-cluster L ----
# Could it be capturing multiple cell type? 
setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/tsom_all_pats/")
dat <- lapply(c(1:4), function(t) {
  fread(paste0("Clustered__Timepoint_T", t, ".csv"))
})
dat <- rbindlist(dat)

dat_meta_l <- dat[dat$TrackSOM_metacluster_lineage_tracking == 'L']

library(Spectre)
dat_meta_l <- run.fitsne(dat_meta_l,
                         use.cols = names(dat_meta_l)[1:38])
make.colour.plot(dat = dat_meta_l,
                 x.axis = "FItSNE_X",
                 y.axis = "FItSNE_Y",
                 col.axis = 'CD3')

markers <- c("CDK1","CD57", "CD274", "CD159c_NKG2C", "CD161", "FOXP3", "CD34", "CD38", "HCMV_pp65", "CD3", "CD278", "CD304", "CD163", "CD314_NKG2D", "CD86", "CD33", "CD14", "CD127", "CD11C", "CD279_PD1", "CD56", "CD19", "CD45RA", "CD69", "CD4", "CD8A", "CD16", "TCRgd", "CD27", "CD20", "Tim3", "CD45RO", "HCMV_IE", "CD25", "Granzyme_B", "CD197_CCR7", "TCRva7.2", "HLA_DR")
func_markers <- c("HCMV_IE", "HCMV_pp65", "CD25", "CD69", "HLA_DR", 'CD38', "CD279_PD1", "Tim3",
                  "CD278", "Granzyme_B", "CD314_NKG2D", "CD86", "CD274")
state_markers <- setdiff(markers, func_markers)
dat_meta_l_stable <- dat_meta_l[, ..state_markers]
dat_meta_l_stable <- melt(dat_meta_l_stable, measure.vars = state_markers)

p <- ggplot(dat_meta_l_stable, aes(x=value)) + 
  # geom_histogram(bins = 100) + 
  geom_density() +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  labs(x = 'Marker expression',
       y = 'Proportion of cells') +
  ggtitle("Expression of cell type markers of cells in meta-cluster L") +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        title = element_text(size=16),
        strip.text = element_text(size=16))
ggsave("meta_l.png",
       height = 12,
       width = 20)

