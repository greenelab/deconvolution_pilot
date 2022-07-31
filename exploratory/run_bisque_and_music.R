# This script is to help determine if deconvolution is possible using the
# pooled single-cell data as reference (i.e. if we can pool future samples).
# We chose to try out Bisque and MuSiC first because a) they're in R and 
# b) they only require cell labels, not creation of a marker gene matrix.
# They also take the exact same inputs, bulk data and single-cell data as
# ExpressionSet objects, so I ran them in the same script.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(BisqueRNA)
  library(MuSiC)
})

source("../config.R")

# Load single cell data
sce <- readRDS("../sce_objects/pooled_clustered_50.rds")

# Load celltypist results
ct <- fread("../celltypist/pooled_predicted_labels.csv")

# Consolidate groups for easier deconvolution
label_table <- fread("../celltypist/simplified_labels.tsv")
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)

# Add cell types to sce object, make changes based on manual annotation
sce$cellType <- ct$Simplified
sce[, sce$cellType == "NK cells" &
    sce$clusters %in% c(2, 6, 10)]$cellType <- "Fibroblasts"
sce[, sce$cellType == "T cells" &
    sce$clusters %in% c(7, 9)]$cellType <- "Epithelial cells"

plotUMAP(sce, colour_by = "cellType")

# Load single cell data into ExpressionSet object
rownames(sce) <- rowData(sce)$ID
sce$SubjectName <- paste(sce$Pool, sce$Barcode, sep = "-")
phenos <- as.data.frame(subset(colData(sce),
                        select = c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription = c("Characters", "Characters"))
pheno <- AnnotatedDataFrame(data = phenos, varMetadata = pheno_metadata)
single_cell <- ExpressionSet(as.matrix(assay(sce)), phenoData = pheno)

# Load bulk data into ExpressionSet object
bulk_matrix <- matrix()
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id,
                         "chunk_ribo/STAR/ReadsPerGene.out.tab", sep = "/")
  bulk_tmp <- fread(bulk_tmp_file)
  bulk_tmp <- bulk_tmp[grep("ENSG", bulk_tmp$V1), ]
  setnames(bulk_tmp, "V2", sample_id)
  if (i == 1) {
    bulk_matrix <- as.matrix(bulk_tmp[, 2])
    rownames(bulk_matrix) <- bulk_tmp$V1
  } else {
    bulk_matrix <- cbind(bulk_matrix, as.matrix(bulk_tmp[, 2]))
  }
}
bulk <- ExpressionSet(bulk_matrix)

# Run bisque
res <- ReferenceBasedDecomposition(bulk.eset = bulk,
                                   sc.eset = single_cell,
                                   use.overlap = FALSE)
# TODO: write this to file once I have a plan for organizing the results
res$bulk.props

# Run music
mus <- music_prop(bulk.eset = bulk,
                  sc.eset = single_cell,
                  clusters = "cellType",
                  samples = "SubjectName")
mus$Est.prop.weighted
mus$Est.prop.allgene

# Load filtered single cell data into ExpressionSet object
fil <- readRDS("../sce_objects/manually_assigned.rds")
rownames(fil) <- rowData(fil)$ID
fil$SubjectName <- paste(fil$Pool, fil$Barcode, sep = "-")
setnames(sce, "cell_type", "cellType")
phenos <- as.data.frame(subset(colData(fil),
                        select = c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription = c("Characters", "Characters"))
pheno <- AnnotatedDataFrame(data = phenos, varMetadata = pheno_metadata)
filter_cell <- ExpressionSet(as.matrix(assay(fil)), phenoData = pheno)

# Run bisque with filtered cell data
res <- ReferenceBasedDecomposition(bulk.eset = bulk,
                                   sc.eset = filter_cell,
                                   use.overlap = FALSE)
res$bulk.props

# Run music with filtered cell data
mus <- music_prop(bulk.eset = bulk,
                  sc.eset = filter_cell,
                  clusters = "cellType",
                  samples = "SubjectName")
mus$Est.prop.weighted
mus$Est.prop.allgene
