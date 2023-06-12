# MuSiC (DOI 10.1038/s41467-018-08023-x)
# https://xuranw.github.io/MuSiC/articles/MuSiC.html

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(MuSiC)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
reference_setting <- snakemake@wildcards[["reference_setting"]]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Get appropriate single cell info and associated output directory
if (is.null(reference_setting)) {
  scefile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
  outdir <- paste(local_data_path, "deconvolution_output", bulk_type, sep = "/")
} else {
  scefile <- paste(local_data_path, "/deconvolution_input/", "labeled_single_cell_profile_", reference_setting, ".rds", sep = "")
  outdir <- paste(local_data_path, "/deconvolution_output/", bulk_type, "_reference_", reference_setting, sep = "")
}

# Load single cell data into ExpressionSet object
# Note: local_data_path is loaded from config.R
sce <- readRDS(scefile)
rownames(sce) <- rowData(sce)$ID
sce$SubjectName <- paste(sce$Pool, sce$Barcode, sep = "-")
phenos <- as.data.frame(subset(colData(sce),
                               select = c("SubjectName", "cellType")))
pheno_metadata <- data.frame(labelDescription = c("Characters", "Characters"))
pheno <- AnnotatedDataFrame(data = phenos, varMetadata = pheno_metadata)
single_cell <- ExpressionSet(as.matrix(assay(sce)), phenoData = pheno)

# Load bulk data into ExpressionSet object
bulk_matrix <- fread(paste(local_data_path, "/deconvolution_input/",
                           "bulk_data_", bulk_type, ".tsv", sep = ""),
                     header = TRUE)
genes <- bulk_matrix$Gene; bulk_matrix$Gene <- NULL
sample_names <- colnames(bulk_matrix)
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes
bulk <- ExpressionSet(bulk_matrix)

# Run music
mus <- music_prop(bulk.eset = bulk,
                  sc.eset = single_cell,
                  clusters = "cellType",
                  samples = "SubjectName")


# Save music object for later perusal
object_file <- paste(outdir, "music_results_full.rds", sep = "/")
saveRDS(mus, file = object_file)

# Format text versions of proportion estimates
nnls <- t(mus$Est.prop.allgene)
nnls <- cbind(rownames(nnls), nnls)
colnames(nnls) <- c("cell_type", sample_names)

music <- t(mus$Est.prop.weighted)
music <- cbind(rownames(music), music)
colnames(music) <- c("cell_type", sample_names)

# Save text versions of proportion estimates
text_file <-  paste(outdir, "music_results.tsv", sep = "/")
write.table(music, file = text_file, sep = "\t", row.names = F, quote = F)

text_file <- paste(outdir, "nnls_results.tsv", sep = "/")
write.table(nnls, file = text_file, sep = "\t", row.names = F, quote = F)
