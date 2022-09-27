# In the BayesPrism tutorial, they specify between "cell type" and "cell state"
# but admit the distinction is murky. They break up their myeloid cells and
# cancer cells, so I will do the same, but also break up the T cells since we
# have so many of them.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(igraph)
  library(dplyr)
})

source("../../config.R")

# Load single cell data
infile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
sce <- readRDS(infile)
sce$unique_barcode <- paste(sce$Pool, sce$Barcode, sep = "-")

# Load celltypist results which have more details as cell state
ct <- fread(paste(local_data_path, "celltypist_output/pooled_predicted_labels.csv", sep = "/"))
ct <- ct[ct$V1 %in% sce$unique_barcode,]
sce$cellState <- ct$majority_voting

# Group macrophages, monocytes, DCs into "myeloid" for cellType
sce$cellType[sce$cellType %in% c("DC","pDC","Monocytes","Macrophages")] <- "Myeloid"

table(sce$cellType)
table(sce$cellState)

# Load in cell labels from genetic demultiplexing
labels_dec <- fread(paste(data_path,"pooled_tumors/12162021/vireo/donor_ids.tsv", sep = "/"))
labels_dec$unique_barcode <- paste("12162021", labels_dec$cell, sep="-")
labels_jan <- fread(paste(data_path,"pooled_tumors/01132022/vireo/donor_ids.tsv", sep = "/"))
labels_jan$unique_barcode <- paste("01132022", labels_jan$cell, sep="-")

# Swap donor ids for the January pool to allow for concatenation
labels_jan$donor_id <- recode(labels_jan$donor_id,
                        "donor0" = "donor4",
                        "donor1" = "donor5",
                        "donor2" = "donor6",
                        "donor3" = "donor7")
labels <- rbind(labels_dec, labels_jan)

# Add donor id to single cell data
labels <- labels[labels$unique_barcode %in% sce$unique_barcode, ]
sce$donor_id <- labels$donor_id

# Switch epithelial cell state to donor info
sce$cellState <- ifelse(sce$cellType=="Epithelial cells", sce$donor_id, sce$cellState)

# Save object
outfile <- paste(local_data_path, "deconvolution_input", "labeled_cell_state_profile.rds", sep = "/")
saveRDS(sce, outfile)