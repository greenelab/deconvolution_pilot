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
})

source("../../config.R")

# Make matrix of pseudobulk with only the miQC-preserved cells
# Note: variable samples is loaded from config.R
filtered_pseudobulk <- matrix(nrow = 36601, ncol = length(samples))
for(i in 1:length(samples)){
  sce <- readRDS(paste("../sce_objects/", samples[i], ".rds", sep=""))
  counts <- as.matrix(assay(sce))
  sums <- rowSums(counts)
  filtered_pseudobulk[, i] <- sums
}
rownames(filtered_pseudobulk) <- names(sums)
colnames(filtered_pseudobulk) <- samples
rm(i, sce, counts, sums); gc()

pseudobulk_path <- paste(local_data_path, "pseudobulk_objects", sep = "/")
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


# Make matrix of pseudobulk but filtered according to Seurat parameters
seurat_pseudobulk <- matrix(nrow = 36601, ncol = length(samples))
for(i in 1:length(samples)){
  sample_id <- samples[i]
  sce <- read10xCounts(paste(data_path, "tumors", sample_id, "Cellranger/outs/filtered_feature_bc_matrix", sep="/"))
  rownames(sce) <- rowData(sce)$Symbol
  mt_genes <- grep("^MT-", rownames(sce),value=T)
  feature_ctrls <- list(mito=mt_genes)
  sce <- addPerCellQC(sce, subsets=feature_ctrls)
  ncol(sce)
  sce <- sce[, sce$subsets_mito_percent<=5]
  sce <- sce[, sce$detected >= 200 & sce$detected <= 2500]
  ncol(sce)
  sce <- logNormCounts(sce)
}
