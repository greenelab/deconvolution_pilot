# XCell (10.1186/s13059-017-1349-1)
# https://github.com/dviraran/xCell


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
bulk_matrix <- fread(bulkfile)
genes <- bulk_matrix$Gene; bulk_matrix$Gene <- NULL
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes

# Run xcell
res <- deconvolute(bulk_matrix, "xcell")

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "xcell_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
