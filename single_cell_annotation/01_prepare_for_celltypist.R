# There's no R implementation of celltypist (nor is there a plan to make one),
# but since the rest of the analysis is so R-centric, I'll convert the sce 
# objects (which have been miQC filtered) to matrices that celltypist can read
# into its python implementation (called run_celltypist.py)

suppressPackageStartupMessages({
  library(data.table)
  library(scater)
  library(scran)
  library(batchelor)
  library(igraph)
})

# This can be run on an individual sample (valid ids 2251, 2267, 2283, 2293, 2380, 2428, 2467, 2497),
# an individual pooled run (valid ids 12162021 and 01132022), or the two pooled runs combined (id should be "pooled")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Sample ID must be provided.", call.=FALSE)
} else {
  sample_id <- args[1]
}

if (sample_id == "pooled") {
  pool1 <- readRDS("sce_objects/12162021.rds")
  pool1$Pool <- "12162021"
  pool2 <- readRDS("sce_objects/01132022.rds")
  pool2$Pool <- "01132022"
  
  sce <- cbind(pool1, pool2)
  
  set.seed(531)
  sce <- runUMAP(sce,
                 BNPARAM = BiocNeighbors::AnnoyParam(),
                 BPPARAM = BiocParallel::MulticoreParam(),
                 min_dist = 0.5,  repulsion_strength = 0.25,
                 spread = 0.7,
                 n_neighbors = 15)
  
  set.seed(531)
  mnn <- fastMNN(sce,  batch = sce$Pool,
                 BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                 BNPARAM = BiocNeighbors::AnnoyParam(),
                 BPPARAM = BiocParallel::MulticoreParam())
  
  reducedDim(sce,  'MNN') <- reducedDim(mnn,  'corrected')
  rm(mnn)
  
  set.seed(531)
  sce <- runUMAP(sce,
                 dimred = 'MNN',  
                 BNPARAM = BiocNeighbors::AnnoyParam(),
                 BPPARAM = BiocParallel::MulticoreParam(),
                 min_dist = 0.5,  repulsion_strength = 0.25,
                 spread = 0.7,
                 n_neighbors = 15)
  
  colnames(sce) <- paste(sce$Pool, sce$Barcode, sep="-")
} else {
  sce <- readRDS(paste("sce_objects/", sample_id, ".rds", sep=""))
  colnames(sce) <- sce$Barcode
}

input_matrix <- as.matrix(assay(sce))
colnames(input_matrix) <- colnames(sce)

input_matrix <- t(input_matrix)

write.table(input_matrix, file=paste("celltypist/", sample_id, "_input_data.csv", sep=""), quote = F, sep=",")