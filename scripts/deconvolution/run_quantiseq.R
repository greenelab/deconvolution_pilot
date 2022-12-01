# quanTIseq (10.1186/s13073-019-0638-6)
# https://icbi.i-med.ac.at/software/quantiseq/doc/

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(immunedeconv)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load salmon results with TPM values
bulkfile <- paste(local_data_path, "/deconvolution_input/normalized_data_", bulk_type, ".tsv", sep = "")
bulk_matrix <- fread(bulkfile, header = TRUE)
genes <- bulk_matrix$Gene; bulk_matrix$Gene <- NULL
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes

# Run quantiseq
res <- deconvolute(bulk_matrix, "quantiseq", tumor = TRUE)

# Consolidate labels based on single-cell categories
labels <- fread(paste(local_data_path,"deconvolution_input",
                      "simplified_labels.tsv", sep = "/"))
setnames(labels, "Original", "cell_type")
res <- inner_join(labels, res)

# Keep granular
res[res$cell_type=="T cell CD4+ (non-regulatory)", ]$Simplified <- "CD4 T cells"
res[res$cell_type=="T cell CD8+", ]$Simplified <- "CD8 T cells"

fractions <- setdiff(colnames(res), c("cell_type","Simplified"))
res <- res %>% group_by(Simplified) %>% summarize_at(vars(fractions), sum)
setnames(res, "Simplified", "cell_type")


# Save quantiseq proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "quantiseq_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
