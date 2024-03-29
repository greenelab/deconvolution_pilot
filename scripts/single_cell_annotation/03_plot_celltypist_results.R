# With cellTypist run, plot the results on a UMAP to see if they look reasonable

suppressPackageStartupMessages({
  library(data.table)
  library(scater)
  library(batchelor)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

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
sce_path <- paste(local_data_path, "sce_objects", sep = "/")
if (sample_id == "pooled") {
  # Load two sets of pools
  pool1 <- readRDS(paste(sce_path, "12162021.rds", sep = "/"))
  pool1$pool <- "12162021"
  pool2 <- readRDS(paste(sce_path, "01132022.rds", sep = "/"))
  pool2$pool <- "01132022"

  sce <- cbind(pool1, pool2)

  # Run batch correction so the two pools will look right on a UMAP
  set.seed(531)
  mnn <- fastMNN(sce, batch = sce$pool,
          BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
          BNPARAM = BiocNeighbors::AnnoyParam(),
          BPPARAM = BiocParallel::MulticoreParam())

  reducedDim(sce, "MNN") <- reducedDim(mnn, "corrected")
  rm(mnn)

  # Run UMAP on batch corrected data
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
  # All other sample_ids will be single units so no need for batch correction.
  # They'll also have had runUMAP generated on their data already.
  sce <- readRDS(paste(sce_path, "/", sample_id, ".rds", sep = ""))
  colnames(sce) <- sce$Barcode
}

# Load celltypist results
celltypist <- fread(paste(local_data_path, "/celltypist_output/", sample_id,
                           "_predicted_labels.csv", sep = ""))
sce$CT_results <- celltypist$majority_voting

# Plot all assignments
png(paste("../../plots/", sample_type, "/", sample_id,
           "_celltypist_full.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "CT_results") + ggtitle(sample_id)
dev.off()

# Consolidate based on hierarchical structure from CellTypist encyclopedia
label_table <- fread(paste(local_data_path,
                           "/celltypist_output/simplified_labels.tsv",
                           sep = ""))
setnames(label_table, "Original", "majority_voting")
celltypist <- left_join(celltypist, label_table)
sce$CT_simplified <- celltypist$Simplified

# Plot simplified assignments
png(paste("../../plots/", sample_type, "/", sample_id,
           "_celltypist_simplified.png", sep = ""), width = 700)
plotUMAP(sce, colour_by = "CT_simplified") + ggtitle(sample_id)
dev.off()
