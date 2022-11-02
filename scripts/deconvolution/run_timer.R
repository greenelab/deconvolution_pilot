# ABIS (10.1007/978-1-0716-0327-7_18)
# http://cistrome.org/TIMER/


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

# Remove the 10 genes with duplicate names
bulk_matrix <- bulk_matrix[which(!duplicated(rownames(bulk_matrix))),]

# Run TIMER
indications <- rep("OV", ncol(bulk_matrix))
res <- deconvolute(bulk_matrix, "timer", indications = indications)

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "timer_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
