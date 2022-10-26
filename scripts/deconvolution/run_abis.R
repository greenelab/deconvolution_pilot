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
bulk_matrix <- fread(bulkfile)
genes <- bulk_matrix$V1; bulk_matrix$V1 <- NULL
bulk_matrix <- as.matrix(bulk_matrix)
rownames(bulk_matrix) <- genes

# Run abis
res <- deconvolute(bulk_matrix, "abis")

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

res <- res %>% group_by(Simplified) %>%
    summarize(`2251`=sum(`2251`),
              `2267`=sum(`2267`),
              `2283`=sum(`2283`),
              `2293`=sum(`2293`),
              `2380`=sum(`2380`),
              `2428`=sum(`2428`),
              `2467`=sum(`2467`),
              `2497`=sum(`2497`))
setnames(res, "Simplified", "cell_type")

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "abis_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
