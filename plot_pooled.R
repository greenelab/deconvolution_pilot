# A first look at the pooled data and the hash demultiplexing done with cellranger multi.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(scater)
  library(miQC)
  library(ggplot2)
})

sample_id <- "12162021"

# Load assignment info for each cell
data_path <- normalizePath(paste("~/Documents/scRNA/sc-cancer-hgsc/data/pooled_tumors", sample_id, sep = "/"))
pool <- fread(paste(data_path, "Cellranger/outs/multi/multiplexing_analysis/assignment_confidence_table.csv", sep = "/"))
table(pool$Assignment)

# Based on this assignment confidence table, the default cutoff for cellranger multi is that it only assigns a cell if
# it has >90% probability of belonging to a given hashtag, <10% probability of being a multiplet, and <10% probability
# of being a blank. Based on a manual look at the assignment table, this seems unnecessarily stringent. Here we test
# what number of cells we get at two different cutoff levels, 85% hashtag and 20% multiplet/blank (moderate), and 90%
# hashtag and 30% multiplet/blank (liberal).
# TODO: determine if the <0.1% requirement for all other hashtags is necessary/helpful. 
moderate <- pool
moderate[moderate$Assignment == "Unassigned" & moderate$`anti-human_Hashtag1` >= .85 & moderate$`anti-human_Hashtag2` <= 0.001 &
           moderate$`anti-human_Hashtag3` <= 0.001 & moderate$`anti-human_Hashtag4` <= 0.001 & moderate$Blanks <= .2 &
           moderate$Multiplet <= .2, ]$Assignment = "anti-human_Hashtag1"
moderate[moderate$Assignment == "Unassigned" & moderate$`anti-human_Hashtag1` <= 0.001 & moderate$`anti-human_Hashtag2` >= .85 &
           moderate$`anti-human_Hashtag3` <= 0.001 & moderate$`anti-human_Hashtag4` <= 0.001 & moderate$Blanks <= .2 &
           moderate$Multiplet <= .2, ]$Assignment = "anti-human_Hashtag2"
moderate[moderate$Assignment == "Unassigned" & moderate$`anti-human_Hashtag1` <= 0.001 & moderate$`anti-human_Hashtag2` <= 0.001 &
           moderate$`anti-human_Hashtag3` >= .85 & moderate$`anti-human_Hashtag4` <= 0.001 & moderate$Blanks <= .2 &
           moderate$Multiplet <= .2, ]$Assignment = "anti-human_Hashtag3"
moderate[moderate$Assignment == "Unassigned" & moderate$`anti-human_Hashtag1` <= 0.001 & moderate$`anti-human_Hashtag2` <= 0.001 &
           moderate$`anti-human_Hashtag3` <= 0.001 & moderate$`anti-human_Hashtag4` >= .85 & moderate$Blanks <= .2 &
           moderate$Multiplet <= .2, ]$Assignment = "anti-human_Hashtag4"
table(moderate$Assignment)

liberal <- pool
liberal[liberal$Assignment == "Unassigned" & liberal$`anti-human_Hashtag1` >= .8 & liberal$`anti-human_Hashtag2` <= 0.001 &
           liberal$`anti-human_Hashtag3` <= 0.001 & liberal$`anti-human_Hashtag4` <= 0.001 & liberal$Blanks <= .3 &
           liberal$Multiplet <= .3, ]$Assignment = "anti-human_Hashtag1"
liberal[liberal$Assignment == "Unassigned" & liberal$`anti-human_Hashtag1` <= 0.001 & liberal$`anti-human_Hashtag2` >= .8 &
           liberal$`anti-human_Hashtag3` <= 0.001 & liberal$`anti-human_Hashtag4` <= 0.001 & liberal$Blanks <= .3 &
           liberal$Multiplet <= .3, ]$Assignment = "anti-human_Hashtag2"
liberal[liberal$Assignment == "Unassigned" & liberal$`anti-human_Hashtag1` <= 0.001 & liberal$`anti-human_Hashtag2` <= 0.001 &
           liberal$`anti-human_Hashtag3` >= .8 & liberal$`anti-human_Hashtag4` <= 0.001 & liberal$Blanks <= .3 &
           liberal$Multiplet <= .3, ]$Assignment = "anti-human_Hashtag3"
liberal[liberal$Assignment == "Unassigned" & liberal$`anti-human_Hashtag1` <= 0.001 & liberal$`anti-human_Hashtag2` <= 0.001 &
           liberal$`anti-human_Hashtag3` <= 0.001 & liberal$`anti-human_Hashtag4` >= .8 & liberal$Blanks <= .3 &
           liberal$Multiplet <= .3, ]$Assignment = "anti-human_Hashtag4"
table(liberal$Assignment)


# Load in total matrix and subset to only the barcodes Cellranger decided were cells
full_matrix <- read10xCounts(paste(data_path, "Cellranger/outs/multi/count/raw_feature_bc_matrix", sep = "/"))
sce <- full_matrix[, full_matrix$Barcode %in% pool$Barcodes]
rm(full_matrix); gc()

# Removing the barcode counts, because they will probably mess up the UMAP
barcodes <- sce[grep("Hashtag", rownames(sce))]
sce <- sce[grep("ENSG", rownames(sce))]

# Prep for miQC
rownames(sce) <- rowData(sce)$Symbol
mt_genes <- grepl("^MT-", rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sce <- addPerCellQC(sce, subsets = feature_ctrls)

# Run miQC
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
               ## unnecessary options,  only used to make a pretty graph
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste("plots/pooled/", sample_id, "_mito.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "subsets_mito_percent")
dev.off()

# Add assignment labels to sce object
assignment <- pool[order(pool$Barcodes), ]
assignment <- subset(assignment, assignment$Barcodes %in% sce$Barcode)$Assignment
assignment[assignment == "Blanks" | assignment == "Multiplet"] <- "Unassigned"
colData(sce) <- cbind(colData(sce), assignment)

# Since one of the samples failed in one of the pools, it has no cells at more stringent thresholds and this
# messes up the indexing for color. Checking for this beforehand each time keeps unassigned cells as gray.
get_colors <- function(assignment) {
  if (length(unique(assignment)) == 5) {
    colors <- scale_color_manual(name = "assignment", values = c("#CA0020", "#0571B0", "#4DAC26", "#E66101", "#999999"))
  } else if (length(unique(assignment)) == 4) {
    colors <- scale_color_manual(name = "assignment", values = c("#0571B0", "#4DAC26", "#E66101", "#999999"))
  }
  colors
}

# Plot UMAP by assignment
png(paste("plots/pooled/", sample_id, "_UMAP_assignment_90.png", sep = ""), width = 700)
colors <- get_colors(sce$assignment); plotUMAP(sce, colour_by = "assignment") + colors
dev.off()

# Check how assignments look on UMAP with other filtering thresholds
mod_assignment <- subset(moderate, moderate$Barcodes %in% sce$Barcode)$Assignment
mod_assignment[mod_assignment == "Blanks" | mod_assignment == "Multiplet"] = "Unassigned"
colData(sce) <- cbind(colData(sce), mod_assignment)
png(paste("plots/pooled/", sample_id, "_UMAP_assignment_85.png", sep = ""), width = 700)
colors <- get_colors(sce$mod_assignment); plotUMAP(sce, colour_by = "mod_assignment") + colors
dev.off()

lib_assignment <- subset(liberal, liberal$Barcodes %in% sce$Barcode)$Assignment
lib_assignment[lib_assignment == "Blanks" | lib_assignment == "Multiplet"] = "Unassigned"
colData(sce) <- cbind(colData(sce), lib_assignment)
png(paste("plots/pooled/", sample_id, "_UMAP_assignment_80.png", sep = ""), width = 700)
colors <- get_colors(sce$lib_assignment); plotUMAP(sce, colour_by = "lib_assignment") + colors
dev.off()


# Check for basic cell type markers (CD45 for immune, EPCAM/PAX8 for epithelial, smooth muscle actin for fibroblasts)
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
