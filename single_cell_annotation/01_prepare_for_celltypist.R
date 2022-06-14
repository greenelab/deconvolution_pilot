# There's no R implementation of celltypist (nor is there a plan to make one),
# but since the rest of the analysis is so R-centric, I'll convert the sce
# objects (which have been miQC filtered) to matrices that celltypist can read
# into its python implementation.

suppressPackageStartupMessages({
  library(data.table)
  library(scater)
  library(scran)
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

# Load sce object
if (sample_id == "pooled") {
  pool1 <- readRDS("../sce_objects/12162021.rds")
  pool1$pool <- "12162021"
  pool2 <- readRDS("../sce_objects/01132022.rds")
  pool2$pool <- "01132022"

  sce <- cbind(pool1, pool2)

  colnames(sce) <- paste(sce$pool, sce$Barcode, sep = "-")
} else {
  sce <- readRDS(paste("../sce_objects/", sample_id, ".rds", sep = ""))
  colnames(sce) <- sce$Barcode
}

# Get count matrix, transpose for celltypist which expects a cell x gene matrix
input_matrix <- as.matrix(assay(sce))
colnames(input_matrix) <- colnames(sce)
input_matrix <- t(input_matrix)

# Write count matrix
write.table(input_matrix, file =
            paste("../celltypist/", sample_id, "_input_data.csv", sep = ""),
            quote = F, sep = ",")
