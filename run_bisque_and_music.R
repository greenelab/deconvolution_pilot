# This script is to help determine if deconvolution is possible using the
# pooled single-cell data as reference (i.e. if we can pool samples in the)
# future. We chose to try out Bisque and MuSiC first because a) they're in R
# and b) they only require cell labels, not creation of a marker gene matrix.
# They also take the exact same inputs, bulk data and single-cell data as
# ExpressionSet objects, so I ran them in the same script.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(BisqueRNA)
  library(MuSiC)
})

source("config.R")

# Load single cell data
sce <- readRDS("sce_objects/pooled_clustered_50.rds")

# Load celltypist results
# TODO: switch to python output
ct <- fread("~/Downloads/celltypist_predictions/predicted_labels.csv")

# Simplify cell types for clearer deconvolution
ct$simplified <- ct$majority_voting
#ct[ct$simplified=="Classical monocytes",]$simplified <- "Monocytes"
#ct[ct$simplified=="pDC",]$simplified <- "DC"

### More simplified
# TODO: switch to table simplification
ct[ct$simplified=="Tem/Trm cytotoxic T cells",]$simplified <- "T cells"
ct[ct$simplified=="Tcm/Naive helper T cells",]$simplified <- "T cells"
ct[ct$simplified=="Regulatory T cells",]$simplified <- "T cells"
ct[ct$simplified=="Classical monocytes",]$simplified <- "Macrophages"
ct[ct$simplified=="Monocytes",]$simplified <- "Macrophages"
ct[ct$simplified=="pDC",]$simplified <- "DC"
ct[ct$simplified=="CD16+ NK cells",]$simplified <- "NK cells"


# Add cell types to sce object, make changes based on manual annotation
sce$cellType <- ct$simplified
sce[, sce$cellType=="NK cells" & sce$clusters %in% c(2, 6, 10)]$cellType <- "Fibroblasts"
sce[, sce$cellType=="Tcm/Naive helper T cells" & sce$clusters %in% c(7, 9)]$cellType <- "Epithelial cells"

plotUMAP(sce, colour_by="cellType")

# Load single cell data into ExpressionSet object
rownames(sce) <- rowData(sce)$ID
sce$SubjectName <- paste(sce$Pool, sce$Barcode, sep = "-")
phenos <- as.data.frame(subset(colData(sce), select=c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription=c("Characters","Characters"))
pheno <- AnnotatedDataFrame(data=phenos, varMetadata = pheno_metadata)
single_cell <- ExpressionSet(as.matrix(assay(sce)), phenoData = pheno)

# Load bulk data into ExpressionSet object
bulk_matrix <- matrix()
samples <- c("2251","2267","2283","2293","2380","2428","2467","2497")
for (i in 1:length(samples)){
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id, "chunk_ribo/STAR/ReadsPerGene.out.tab", sep="/")
  bulk_tmp <- fread(bulk_tmp_file)
  bulk_tmp <- bulk_tmp[grep("ENSG", bulk_tmp$V1), ]
  setnames(bulk_tmp, "V2",sample_id)
  if (i == 1){
    bulk_matrix <- as.matrix(bulk_tmp[, 2])
    rownames(bulk_matrix) <- bulk_tmp$V1
  } else {
    bulk_matrix <- cbind(bulk_matrix, as.matrix(bulk_tmp[, 2]))
  }
}
bulk <- ExpressionSet(bulk_matrix)

# Run bisque
res <- ReferenceBasedDecomposition(bulk.eset = bulk, sc.eset = single_cell, use.overlap = F)
res$bulk.props

# Run music
mus <- music_prop(bulk.eset = bulk, sc.eset = single_cell, clusters="cellType", samples="SubjectName")
mus$Est.prop.weighted
mus$Est.prop.allgene