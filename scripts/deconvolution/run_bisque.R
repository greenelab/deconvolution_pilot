# bisque (DOI 10.1038/s41467-020-15816-6)
# https://github.com/cozygene/bisque

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(BisqueRNA)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load single cell data into ExpressionSet object
# Note: local_data_path is loaded from config.R
sce <- readRDS(paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/"))
rownames(sce) <- rowData(sce)$ID
sce$SubjectName <- paste(sce$Pool, sce$Barcode, sep = "-")
phenos <- as.data.frame(subset(colData(sce),
                               select = c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription = c("Characters", "Characters"))
pheno <- AnnotatedDataFrame(data = phenos, varMetadata = pheno_metadata)
single_cell <- ExpressionSet(as.matrix(assay(sce)), phenoData = pheno)

# Load bulk data into ExpressionSet object
bulk_matrix <- fread(paste(local_data_path, "/deconvolution_input/",
                           "bulk_data_", bulk_type, ".tsv", sep = ""))
genes <- bulk_matrix$V1; bulk_matrix$V1 <- NULL
sample_names <- colnames(bulk_matrix)
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes
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
colnames(tmp) <- c("cell_type",sample_names)

# Save proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "bisque_results.tsv", sep = "/")
write.table(tmp, file = text_file, sep = "\t", row.names = F, quote=F)
