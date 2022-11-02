# Consensus_TME (10.1158/0008-5472.CAN-18-3560)
# https://github.com/cansysbio/ConsensusTME

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

# Consensus_tme asks for a list of TCGA cancer types that each sample belongs to
indications <- rep("ov", ncol(bulk_matrix))

# Run consensus_tme
res <- deconvolute_consensus_tme(bulk_matrix, indications)
sample_names <- colnames(res)
res <- cbind(rownames(res), res)
colnames(res) <- c("cell_type", sample_names)

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "consensus_tme_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)