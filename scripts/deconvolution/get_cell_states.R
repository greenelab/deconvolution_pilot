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
  library(yaml)
})

reference_setting <- snakemake@wildcards[["reference_setting"]]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load single cell data
if (is.null(reference_setting)){
  infile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
} else{
  infile <- paste(local_data_path, "/deconvolution_input/", "labeled_single_cell_profile_", reference_setting, ".rds", sep = "")
}
sce <- readRDS(infile)
sce$unique_barcode <- paste(sce$Pool, sce$Barcode, sep = "-")

# Load celltypist results which have more details as cell state
ct <- fread(paste(local_data_path, "celltypist_output/pooled_predicted_labels.csv", sep = "/"))
ct <- ct[ct$V1 %in% sce$unique_barcode,]
sce$cellState <- ct$majority_voting

table(sce$cellType)
table(sce$cellState)

# Load in cell labels from genetic referenceing
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
sce$cellState <- ifelse(sce$cellType == "Epithelial cells", sce$donor_id, sce$cellState)

# Save object
if (is.null(reference_setting)){
  outfile <- paste(local_data_path, "deconvolution_input", "labeled_cell_state_profile.rds", sep = "/")
} else {
  outfile <- paste(local_data_path, "/deconvolution_input/", "labeled_cell_state_profile_", reference_setting, ".rds", sep = "")
}
saveRDS(sce, outfile)
