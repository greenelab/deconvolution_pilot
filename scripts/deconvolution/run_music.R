# MuSiC (DOI 10.1038/s41467-018-08023-x)
# https://xuranw.github.io/MuSiC/articles/MuSiC.html

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(MuSiC)
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

# Run music
mus <- music_prop(bulk.eset = bulk,
                  sc.eset = single_cell,
                  clusters = "cellType",
                  samples = "SubjectName")


# Save music object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "music_results_full.rds", sep = "/")
saveRDS(mus, file = object_file)

# Save text versions of proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "music_results.tsv", sep = "/")
write.table(mus$Est.prop.weighted, file = text_file, sep = "\t")

text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "nnls_results.tsv", sep = "/")
write.table(mus$Est.prop.allgene, file = text_file, sep = "\t")