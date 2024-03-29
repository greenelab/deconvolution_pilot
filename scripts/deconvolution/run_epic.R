# EPIC (10.1007/978-1-0716-0327-7_17)
# https://github.com/GfellerLab/EPIC

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(EPIC)
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

# Run epic with default reference cells
out <- EPIC(bulk = bulk_matrix)

# Save epic object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "epic_results_full.rds", sep = "/")
saveRDS(out, file = object_file)

# Reformat text version of proportion estimates
tmp <- as.data.frame(t(out$cellFractions))
tmp <- cbind(rownames(tmp), tmp)
colnames(tmp) <- c("cell_type", colnames(bulk_matrix))

# Rename labels based on single-cell categories
# Based on a line in the Racle and Gfeller paper, we're going to
# list the "otherCells" type as epithelial. TODO: check with Casey
labels <- fread(paste(local_data_path, "deconvolution_input",
                      "simplified_labels.tsv", sep = "/"))
setnames(labels, "Original", "cell_type")
tmp <- inner_join(labels, tmp)

fractions <- setdiff(colnames(tmp), c("cell_type","Simplified"))
tmp <- tmp %>% group_by(Simplified) %>% summarize_at(vars(fractions), sum)
setnames(tmp, "Simplified", "cell_type")

# Save proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "epic_results.tsv", sep = "/")
write.table(tmp, file = text_file, sep = "\t", row.names = F, quote = F)
