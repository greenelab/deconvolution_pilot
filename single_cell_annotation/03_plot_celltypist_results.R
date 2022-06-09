# With cellTypist run, plot the results on a UMAP to see if they look reasonable

suppressPackageStartupMessages({
 library(data.table)
 library(scater)
 library(batchelor)
 library(dplyr)
})

# This can be run on an individual sample (valid ids 2251, 2267, 2283, 2293,
# 2380, 2428, 2467, 2497), an individual pooled run (valid ids 12162021 and
# 01132022), or the two pooled runs combined (valid id pooled).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
 stop("Sample ID must be provided.", call. = FALSE)
} else {
 sample_id <- args[1]
}

# Determine which output directory to put plots in
if (sample_id %in% c("12162021", "01132022", "pooled")) {
 sample_type <- "pooled"
} else {
 sample_type <- "single_cell"
}

# Load sce object, if using all pooled data run batchelor on the
# combined samples to correct for batch effects
if (sample_id == "pooled") {
 pool1 <- readRDS("../sce_objects/12162021.rds")
 pool1$pool <- "12162021"
 pool2 <- readRDS("../sce_objects/01132022.rds")
 pool2$pool <- "01132022"

 sce <- cbind(pool1, pool2)

 set.seed(531)
 mnn <- fastMNN(sce, batch = sce$pool,
         BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
         BNPARAM = BiocNeighbors::AnnoyParam(),
         BPPARAM = BiocParallel::MulticoreParam())

 reducedDim(sce, "MNN") <- reducedDim(mnn, "corrected")
 rm(mnn)

 set.seed(531)
 sce <- runUMAP(sce,
         dimred = "MNN",
         BNPARAM = BiocNeighbors::AnnoyParam(),
         BPPARAM = BiocParallel::MulticoreParam(),
         min_dist = 0.5, repulsion_strength = 0.25,
         spread = 0.7,
         n_neighbors = 15)

 colnames(sce) <- paste(sce$pool, sce$Barcode, sep = "-")
} else {
 sce <- readRDS(paste("../sce_objects/", sample_id, ".rds", sep = ""))
 colnames(sce) <- sce$Barcode
}

# Load celltypist results
celltypist <- fread(paste("../celltypist/", sample_id,
                          "_predicted_labels.csv", sep = ""))
sce$CT_results <- celltypist$majority_voting

# Plot all assignments
png(paste("../plots/", sample_type, "/", sample_id,
          "_cellassign_full.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "CT_results") + ggtitle(sample_id)
dev.off()

# Consolidate based on hierarchical structure from CellTypist encyclopedia
label_table <- fread("../celltypist/simplified_labels.tsv")
setnames(label_table, "Original", "majority_voting")
celltypist <- left_join(celltypist, label_table)
sce$CT_simplified <- celltypist$Simplified

# Plot simplified assignments
png(paste("../plots/", sample_type, "/", sample_id,
          "_cellassign_simplified.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "CT_simplified") + ggtitle(sample_id)
dev.off()
