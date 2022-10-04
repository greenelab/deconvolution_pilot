# Most deconvolution methods that employ single-cell data require a reference
# profile with each cell labeled by cell type. This script generates that,
# using the intersection of cellTypist and unsupervised clustering results
# to determine cell labels. Note: this is using data from the pooled and
# demultiplexed runs, meaing that we can use the individually sequenced single
# cell runs to generate pseudobulk files to evaluate without testing on our
# training data.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
})

source("../../config.R")

# Load single cell data
sce <- readRDS(paste(local_data_path, 
                     "sce_objects/pooled_clustered_50.rds", sep = "/"))

# Load celltypist results
ct <- fread(paste(local_data_path,
                  "celltypist_output/pooled_predicted_labels.csv", sep = "/"))

# Load in cell labels from genetic demultiplexing
labels_dec <- fread(paste(data_path,"pooled_tumors/12162021/vireo/donor_ids.tsv", sep = "/"))
labels_dec$unique_barcode <- paste("12162021", labels_dec$cell, sep="-")
labels_jan <- fread(paste(data_path,"pooled_tumors/01132022/vireo/donor_ids.tsv", sep = "/"))
labels_jan$unique_barcode <- paste("01132022", labels_jan$cell, sep="-")
labels <- rbind(labels_dec, labels_jan)

# Remove doublets and unassigned cells
sce$unique_barcode <- paste(sce$Pool, sce$Barcode, sep = "-")
labels <- labels[labels$unique_barcode %in% sce$unique_barcode, ]
sce$donor_id <- labels$donor_id
sce <- sce[, sce$donor_id %in% c("donor0", "donor1", "donor2", "donor3")]
ct <- ct[ct$V1 %in% sce$unique_barcode, ]

# Consolidate groups for easier deconvolution
label_table <- fread(paste(local_data_path,
                     "celltypist_output/simplified_labels.tsv", sep = "/"))
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)

# Get overlap of cellTypist and manual annotation of unsupervised clustering
ct$keep <- FALSE
ct$clusters <- sce$clusters
ct[ct$clusters %in% c(1, 3, 12) & ct$Simplified %in% c("T cells", "NK cells", "ILC", "Mast cells"),]$keep <- TRUE
ct[ct$clusters %in% c(2, 6, 10) & ct$Simplified=="Fibroblasts", ]$keep <- TRUE
ct[ct$clusters==4 & ct$Simplified=="Endothelial cells",]$keep <- TRUE
ct[ct$clusters==5 & ct$Simplified %in% c("Monocytes", "Macrophages", "DC"),]$keep <- TRUE
ct[ct$clusters %in% c(7, 9) & ct$Simplified=="Epithelial cells",]$keep <- TRUE
ct[ct$clusters==8 & ct$Simplified %in% c("pDC", "B cells"),]$keep <- TRUE
ct[ct$clusters==11 & ct$Simplified=="Plasma cells",]$keep <- TRUE

# Keep cells within overlap
sce$cellType <- ct$Simplified
sce <- sce[, ct$keep]

# Save object
outfile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
saveRDS(sce, file = outfile)
