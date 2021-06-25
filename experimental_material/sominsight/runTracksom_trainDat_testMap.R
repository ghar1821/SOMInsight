# TrackSOM is ran on training data.
# Testing data is mapped onto FlowSOM's SOM nodes.
# Based on the node membership, determine the cell's meta clusters
# Import libraries
library(TrackSOM)
library(Spectre)
library(tools)

data_dir <- "M:/givanna/cmv/data_complete_timepoints_only/"
data.files <- list.files(data_dir, ".csv")
dat_per_tp <- lapply(data.files, function(f) fread(paste0(data_dir, f)))

markers <- c("CDK1","CD57", "CD274", "CD159c_NKG2C", "CD161", "FOXP3", "CD34", "CD38", "HCMV_pp65", "CD3", "CD278", "CD304", "CD163", "CD314_NKG2D", "CD86", "CD33", "CD14", "CD127", "CD11C", "CD279_PD1", "CD56", "CD19", "CD45RA", "CD69", "CD4", "CD8A", "CD16", "TCRgd", "CD27", "CD20", "Tim3", "CD45RO", "HCMV_IE", "CD25", "Granzyme_B", "CD197_CCR7", "TCRva7.2", "HLA_DR")

# We divide the data into folds. So we need to identify which fold to use to extract
# features from
fold_files_dir <- "M:/givanna/cmv/loocv_folds/"
setwd(fold_files_dir)
fold_files <- sapply(list.files(), function(x) paste0(fold_files_dir, x))

save_dir <- "M:/givanna/cmv/tracksom_10x10_20_noMerge_completeTP/"
dir.create(save_dir, recursive = TRUE)
# Setup and run diffyct per fold ----
for (fold_file in fold_files) {
    # fold_file <- fold_files[1]
    fold_split <- fread(fold_file)
    fold_train <- fold_split[fold_split$set == 'train',]
    
    fold_id <- file_path_sans_ext(basename(fold_file))
    save_dir_fold <- paste0(save_dir, fold_id)
    dir.create(save_dir_fold)
    setwd(save_dir_fold)
    
    message(paste("Processing", fold_id))
    
    dat_train <- lapply(dat_per_tp, function(d) {
        d_subset <- d[d$PatientName %in% fold_train$PatientName, ]
        return(d_subset)
    })
    
    # try small subset just for testing
    # dat_train <- lapply(dat_train, function(d) head(d, 10000))
    
    message("Running TrackSOM")
    tracksom.result <- TrackSOM(inputFiles = dat_train,
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
    dat_train <- rbindlist(dat_train)
    dat_train <- ConcatenateClusteringDetails(tracksom.result = tracksom.result,
                                             dat = dat_train,
                                             timepoint.col = 'Timepoint',
                                             timepoints = c("T1", "T2", "T3", "T4"))
    message("Saving training data")
    Spectre::write.files(dat_train, file.prefix = 'Clustered_trainDat', divide.by = 'Timepoint')
    
    message("Mapping test data to FlowSOM code")
    # try to assign cluster to the test data now.
    fold_test <- fold_split[fold_split$set == 'test',]
    
    # map test data
    dat_test <- lapply(dat_per_tp, function(d) {
        d_subset <- d[d$PatientName %in% fold_test$PatientName, ]
        return(d_subset)
    })
    
    dat_test_ff <- lapply(dat_test, function(d) {
        input_dat <- d[, markers, with = FALSE]
        
        # convert to flowFrame. Pain in the butt
        metadata <- data.frame(name=dimnames(input_dat)[[2]], desc=paste('column',dimnames(input_dat)[[2]],'from dataset'))
        input_dat_ff <- new("flowFrame",
                            exprs=as.matrix(input_dat), # in order to create a flow frame, data needs to be read as matrix
                            parameters=Biobase::AnnotatedDataFrame(metadata))
        return(input_dat_ff)
    })
    dat_test_flowset <- as(dat_test_ff, "flowSet")
    
    dat_test <- rbindlist(dat_test)
    
    fsom_test <- FlowSOM::NewData(tracksom.result$FlowSOM, dat_test_flowset)
    
    message("Extracting test data's code and meta cluster's ID")
    fsom_test_codes <- data.table(TrackSOM_cluster=fsom_test$map$mapping[, 1])
    fsom_test_codes$Timepoint <- dat_test$Timepoint
    fsom_test_codes$row_idx <- c(1:nrow(fsom_test_codes))
    
    # Get meta for each code from training dat
    dat_train_metas <- unique(dat_train[, c("TrackSOM_cluster", 
                                            "TrackSOM_metacluster", 
                                            "TrackSOM_metacluster_lineage_tracking", 
                                            "Timepoint")])
    # Just do a merge to get the meta clusters
    fsom_test_codes <- merge(fsom_test_codes, dat_train_metas, all.x = TRUE, all.y = FALSE)
    fsom_test_codes <- fsom_test_codes[order(row_idx)]
    
    # Remove time point and row_idx
    fsom_test_codes[, c("Timepoint", "row_idx") := NULL]
    
    # Assign the metas
    to_assign <- c("TrackSOM_cluster", 
                   "TrackSOM_metacluster", 
                   "TrackSOM_metacluster_lineage_tracking")
    for (t in to_assign) {
        dat_test[[t]] <- fsom_test_codes[[t]]
    }
    
    message("Writing test data out")
    Spectre::write.files(dat_test, file.prefix = "Clustered_testDat", divide.by = 'Timepoint')
}




