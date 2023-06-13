# CibersortX uses single-cell data where the labels are embedded into
# the matrix as colnames. I'll format it here.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(yaml)
})

reference_setting <- snakemake@wildcards[["reference_setting"]]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load in single cell data
if (is.null(reference_setting)) {
  infile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
} else{ 
  infile <- paste(local_data_path, "/deconvolution_input/labeled_single_cell_profile_", reference_setting, ".rds", sep = "")
}
sce <- readRDS(infile)
rownames(sce) <- rowData(sce)$ID

counts <- as.matrix(assay(sce))
new_counts <- cbind(rownames(counts),counts)
colnames(new_counts) <- c("Gene",sce$cellType)

if (is.null(reference_setting)) {
  outfile <- paste(local_data_path, "deconvolution_input", "cibersortx_single_cell_profile.tsv", sep = "/")
} else {
  outfile <- paste(local_data_path, "/deconvolution_input/cibersortx_single_cell_profile_", reference_setting, ".tsv", sep = "")
}
write.table(new_counts, file = outfile, quote=F, row.names=F, sep="\t")
