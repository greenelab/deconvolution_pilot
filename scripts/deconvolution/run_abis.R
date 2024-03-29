# ABIS (10.1016/j.celrep.2019.01.041)
# https://github.com/giannimonaco/ABIS


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

# Run abis
res <- deconvolute(bulk_matrix, "abis")

print(res)

# The abis authors say that small negative values (>-5%) are
# expected as a result of noise and can be set to 0
res <- as.data.frame(res)
if(min(res[,-1]) < -5){
  stop("Error, found a cell type with proportion of -5% or less.
       This is a sign of strong technical or biological variability
       and this cell type should be excluded from analysis.")
}
res[res < 0] <- 0

# Switch to percentage values to allow direct comparison to other proportions
res <- cbind(cell_type=res$cell_type, res[,-1]/100)
  
# Consolidate labels based on single-cell categories
labels <- fread(paste(local_data_path,"deconvolution_input",
                      "simplified_labels.tsv", sep = "/"))
setnames(labels, "Original", "cell_type")
res <- inner_join(labels, res)

# Keep it granular
res[res$cell_type=="T cell CD8+ memory" |
    res$cell_type=="T cell CD8+ naive", ]$Simplified <- "CD8 T cells"
res[res$cell_type=="T cell CD4+ memory" |
    res$cell_type=="T cell CD4+ naive", ]$Simplified <- "CD4 T cells"

fractions <- setdiff(colnames(res), c("cell_type","Simplified"))
res <- res %>% group_by(Simplified) %>% summarize_at(vars(fractions), sum)
setnames(res, "Simplified", "cell_type")

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "abis_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
