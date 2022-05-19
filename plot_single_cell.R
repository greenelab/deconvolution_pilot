# First look at the regular (non-pooled) single cell data from our pilot samples.

setwd("~/Documents/scRNA/sc-cancer-hgsc/data")

suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(scater)
  library(miQC)
  library(ggplot2)
})

sample_id <- "2497"

# Load data
data_path <- normalizePath("~/Documents/scRNA/sc-cancer-hgsc/data/tumors/")
sce <- read10xCounts(paste(data_path, sample_id, "Cellranger/outs/filtered_feature_bc_matrix", sep = "/"))
rownames(sce) <- rowData(sce)$Symbol

# Run miQC
mt_genes <- grepl("MT-", rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sce <- addPerCellQC(sce, subsets = feature_ctrls)
model <- mixtureModel(sce)
p <- plotFiltering(sce, model)

before_cells <- ncol(sce)
sce <- filterCells(sce, model)
after_cells <- ncol(sce)

# Plot miQC results
png(paste("plots/single_cell/", sample_id, "_miQC.png", sep = ""), width = 700)
p + ggtitle(paste(sample_id, "(Keeping", after_cells, "out of", before_cells, "cells)"))
dev.off()

sce <- logNormCounts(sce)

# Make UMAP
set.seed(317)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               ## unnecessary options,  only used to make a pretty graph
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

# Plot a few basic markers. Some of these are markers of cell type (CD45 for immune,
# ACTA2 and VIM for fibroblasts, EPCAM and PAX8 for epithelial), some of them are
# markers of proliferation (MKI67, AURKA, and MYC).
png(paste("plots/single_cell/", sample_id, "_UMAP_PAX8.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "PAX8") + ggtitle(paste(sample_id, "PAX8"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_CD45.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "PTPRC") + ggtitle(paste(sample_id, "CD45"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_EPCAM.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "EPCAM") + ggtitle(paste(sample_id, "EPCAM"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_ACTA2.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "ACTA2") + ggtitle(paste(sample_id, "smooth muscle actin"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_VIM.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "VIM") + ggtitle(paste(sample_id, "VIM"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_MKI67.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "MKI67") + ggtitle(paste(sample_id, "KI67"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_AURKA.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "AURKA") + ggtitle(paste(sample_id, "Aurora kinase"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_UMAP_MYC.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "MYC") + ggtitle(paste(sample_id, "MYC"))
dev.off()

png(paste("plots/single_cell/", sample_id, "_mito.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "subsets_mito_percent") + ggtitle(paste(sample_id, "% Mito"))
dev.off()
