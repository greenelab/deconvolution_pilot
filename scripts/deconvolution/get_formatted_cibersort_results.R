suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

infile <- paste(local_data_path, "deconvolution_output", bulk_type, "CIBERSORTx_Results.txt", sep = "/")
res <- fread(infile)

res$`P-value` <- NULL
res$Correlation <- NULL
res$RMSE <- NULL

sample_names <- res$Mixture
res$Mixture <- NULL

res <- t(res)
res <- cbind(rownames(res),res)
colnames(res) <- c("cell_type", sample_names)

outfile <- paste(local_data_path, "deconvolution_output", bulk_type, "cibersortx_results.tsv", sep = "/")
write.table(res, outfile, sep="\t", quote=F, row.names = F)
