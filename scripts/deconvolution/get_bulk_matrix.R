# Since we're running CIBERSORTx on the console instead of inside a script,
# we need the bulk data matrix saved to file. TODO: decide if I should go
# back and have the other deconvolution methods load in this file.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['truebulktype']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load bulk data
bulk_matrix <- matrix()
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id, bulk_type,
                         "STAR/ReadsPerGene.out.tab", sep = "/")
  bulk_tmp <- fread(bulk_tmp_file)
  bulk_tmp <- bulk_tmp[grep("ENSG", bulk_tmp$V1), ]
  setnames(bulk_tmp, "V2", sample_id)
  if (i == 1) {
    bulk_matrix <- as.matrix(bulk_tmp[, 2])
    rownames(bulk_matrix) <- bulk_tmp$V1
  } else {
    bulk_matrix <- cbind(bulk_matrix, as.matrix(bulk_tmp[, 2]))
  }
}

outfile <- paste(local_data_path, "/deconvolution_input/bulk_data_", bulk_type, ".tsv", sep = "")
write.table(bulk_matrix, outfile, quote=F, sep="\t")
