# CibersortX uses single-cell data where the labels are embedded into
# the matrix as colnames. I'll format it here.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(yaml)
})

demultiplex_setting <- snakemake@wildcards[["demultiplex_setting"]]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load in single cell data
if (is.null(demultiplex_setting)) {
  infile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
  print("got in here okay")
} else{ 
  infile <- paste(local_data_path, "/deconvolution_input/labeled_single_cell_profile_", demultiplex_setting, ".rds", sep = "")
  print("not null")
}
sce <- readRDS(infile)
rownames(sce) <- rowData(sce)$ID

counts <- as.matrix(assay(sce))
print("and here")
new_counts <- cbind(rownames(counts),counts)
print("All good")
colnames(new_counts) <- c("Gene",sce$cellType)
print("Yep yep yep")

if (is.null(demultiplex_setting)) {
  outfile <- paste(local_data_path, "deconvolution_input", "cibersortx_single_cell_profile.tsv", sep = "/")
  print("outfile for null")
} else {
  outfile <- paste(local_data_path, "/deconvolution_input/cibersortx_single_cell_profile_", demultiplex_setting, ".tsv", sep = "")
  print("outfile for not null")
}
write.table(new_counts, file = outfile, quote=F, row.names=F, sep="\t")
