# bisque (DOI 10.1038/s41467-020-15816-6)
# https://github.com/cozygene/bisque

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(BisqueRNA)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
source("../../config.R")

# Load single cell data into ExpressionSet object
sce <- readRDS(paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/"))
rownames(sce) <- rowData(sce)$ID
sce$SubjectName <- paste(sce$Pool, sce$Barcode, sep = "-")
phenos <- as.data.frame(subset(colData(sce),
                               select = c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription = c("Characters", "Characters"))
pheno <- AnnotatedDataFrame(data = phenos, varMetadata = pheno_metadata)
single_cell <- ExpressionSet(as.matrix(assay(sce)), phenoData = pheno)

# Load bulk data into ExpressionSet object
# Note: data_path is loaded from config.R
bulk_matrix <- matrix()
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id, bulk_type,
                         "STAR/ReadsPerGene.out.tab", sep = "/")
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

# Save bisque object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "bisque_results_full.rds", sep = "/")
saveRDS(res, file = object_file)

# Format text version of proportion estimates
tmp <- cbind(rownames(res$bulk.props), res$bulk.props)
colnames(tmp) <- c("cell_type",samples)

# Save proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "bisque_results.tsv", sep = "/")
write.table(tmp, file = text_file, sep = "\t", row.names = F, quote=F)
