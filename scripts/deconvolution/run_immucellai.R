# ImmuCellAI (10.1002/advs.201902880)
# https://github.com/lydiaMyr/ImmuCellAI

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ImmuCellAI)
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
genes <- bulk_matrix$V1; bulk_matrix$V1 <- NULL
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes

# group_tag is 0 because we're not trying to compare between different groups
# response_tag is 0 because we don't need to predict response to ICB therapy
res <- ImmuCellAI_new(bulk_matrix, "rnaseq", group_tag = 0, response_tag = 0)

# Save full object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "immucellai_results_full.rds", sep = "/")
saveRDS(res, file = object_file)

# Reformat output data to match style of other deconvolution methods
out <- t(res$Sample_abundance)
out <- cbind(rownames(out),out)
colnames(out) <- c("cell_type", colnames(bulk_matrix))

outfile <- paste(local_data_path, "deconvolution_output", bulk_type,
                 "immucellai_results.tsv", sep = "/")
write.table(out, outfile, sep = "\t", quote=F, row.names = F)
