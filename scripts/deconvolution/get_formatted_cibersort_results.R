# Since cibersortx isn't run within a script, it needs a separate script to
# rename and format to match the other methods' results.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
})

demultiplex_setting <- snakemake@wildcards[["demultiplex_setting"]]
bulk_type <- snakemake@wildcards[["bulk_type"]]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

if (is.null(demultiplex_setting)) {
  outdir <- paste(local_data_path, "deconvolution_output", bulk_type, sep = "/")
} else {
  outdir <- paste(local_data_path, "/deconvolution_output/", bulk_type, "_demultiplex_default", sep = "")
}

infile <- paste(outdir, "CIBERSORTx_Results.txt", sep = "/")
res <- fread(infile)

res$`P-value` <- NULL
res$Correlation <- NULL
res$RMSE <- NULL

sample_names <- res$Mixture
res$Mixture <- NULL

res <- t(res)
res <- cbind(rownames(res), res)
colnames(res) <- c("cell_type", sample_names)

outfile <- paste(outdir, "cibersortx_results.tsv", sep = "/")
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)