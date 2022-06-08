# A first look at the pooled data from our pilot samples. We will try separating these out by
# donor using both hashing demultiplexing and genetic demultiplexing.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(scater)
  library(miQC)
  library(ggplot2)
})

source("config.R")

# We have two sets of pooled samples, named after the date they were run (12162021 and 01132022).
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Sample ID must be provided.", call.=FALSE)
} else {
  sample_id <- args[1]
}

# Load hashing demultiplexing assignment info to get barcodes to keep
hashing <- fread(paste(data_path, "pooled_tumors", sample_id,
		       "Cellranger/outs/multi/multiplexing_analysis/assignment_confidence_table.csv", sep = "/"))

# Load in total matrix and subset to only the barcodes Cellranger decided were cells
full_matrix <- read10xCounts(paste(data_path, "pooled_tumors", sample_id,
				   "Cellranger/outs/multi/count/raw_feature_bc_matrix", sep = "/"))
sce <- full_matrix[, full_matrix$Barcode %in% hashing$Barcodes]
rm(full_matrix); gc()

# Removing the barcode counts, they will mess with miQC results
barcodes <- sce[grep("Hashtag", rownames(sce))]
sce <- sce[grep("ENSG", rownames(sce))]

# Prep for miQC
rownames(sce) <- rowData(sce)$Symbol
mt_genes <- grepl("^MT-", rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sce <- addPerCellQC(sce, subsets = feature_ctrls)

# Run miQC
set.seed(516)
model <- mixtureModel(sce)
plotModel(sce, model)
p <- plotFiltering(sce, model)

before_cells <- ncol(sce)
sce <- filterCells(sce, model)
after_cells <- ncol(sce)

# Plot miQC results
png(paste("plots/pooled/", sample_id, "_miQC.png", sep = ""), width = 700)
p + ggtitle(paste(sample_id, "(Keeping", after_cells, "out of", before_cells, "cells)"))
dev.off()

# UMAP
set.seed(516)
sce <- logNormCounts(sce)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste("plots/pooled/", sample_id, "_mito.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "subsets_mito_percent")
dev.off()

# Check for basic cell type markers (CD45 (aka PTPRC) for immune, EPCAM/PAX8 for
# epithelial, smooth muscle actin for fibroblasts)
png(paste("plots/pooled/", sample_id, "_UMAP_CD45.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "PTPRC") + ggtitle(paste(sample_id, "CD45"))
dev.off()

png(paste("plots/pooled/", sample_id, "_UMAP_EPCAM.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "EPCAM") + ggtitle(paste(sample_id, "EPCAM"))
dev.off()

png(paste("plots/pooled/", sample_id, "_UMAP_PAX8.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "PAX8") + ggtitle(paste(sample_id, "PAX8"))
dev.off()

png(paste("plots/pooled/", sample_id, "_UMAP_ACTA2.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "ACTA2") + ggtitle(paste(sample_id, "smooth muscle actin"))
dev.off()

# Save filtered and UMAP'ed SingleCellExperiment object for future use
saveRDS(sce, file=paste("sce_objects/", sample_id, ".rds", sep = ""))
