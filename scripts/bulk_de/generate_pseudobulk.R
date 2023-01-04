# We want to use pseudo-bulk data (single cell data where the cells are all
# pooled together to approximate a bulk RNA-seq sample) several ways in this
# study, including running deconvolution on it as a control, and running
# differential expression with it and the true bulk RNA-seq data. We'll make
# two pseudobulk files, one with all cells returned by cellRanger and one
# with only cells that passed miQC filtering.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scater)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples


# Make matrix of pseudobulk with only the miQC-preserved cells
# Note: variable samples is loaded from config.R
filtered_pseudobulk <- matrix(nrow = 36601, ncol = length(samples))
for(i in 1:length(samples)){
  sce <- readRDS(paste(local_data_path, "/sce_objects/", samples[i], ".rds", sep=""))
  counts <- as.matrix(assay(sce))
  sums <- rowSums(counts)
  filtered_pseudobulk[, i] <- sums
}
rownames(filtered_pseudobulk) <- names(sums)
colnames(filtered_pseudobulk) <- samples
rm(i, sce, counts, sums); gc()

pseudobulk_path <- paste(local_data_path, "sce_objects", sep = "/")
saveRDS(filtered_pseudobulk, file = paste(pseudobulk_path, "filtered_pseudobulk.rds", sep = "/"))

# Make matrix of pseudobulk with all cells
# Note: variable data_path is loaded from config.R
full_pseudobulk <- matrix(nrow = 36601, ncol = length(samples))
for(i in 1:length(samples)){
  sample_id <- samples[i]
  sce <- read10xCounts(paste(data_path, "tumors", sample_id, "Cellranger/outs/filtered_feature_bc_matrix", sep = "/"))
  counts <- as.matrix(assay(sce))
  sums <- rowSums(counts)
  full_pseudobulk[, i] <- sums
}
rownames(full_pseudobulk) <- rowData(sce)$Symbol
colnames(full_pseudobulk) <- samples
rm(i, sce, counts, sums); gc()

saveRDS(full_pseudobulk, file = paste(pseudobulk_path, "full_pseudobulk.rds", sep = "/"))
