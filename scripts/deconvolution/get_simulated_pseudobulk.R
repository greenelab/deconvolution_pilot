# We're going to try running all of our deconvolution methods on pseudobulk
# data, as the more expected "gold standard" way of assessing deconvolution
# accuracy, with the caveat that we see some gene length bias we're still
# not totally sure how to account for. SimBu has developed a package that
# makes it easy to generate pseudobulk data from single-cell data, with several
# different possible scenarios. We'll try several here.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(SimBu)
  library(yaml)
})

# There's a bug in the SimBu code which causes mislabeled cell types in one
# use case. I submitted it as an issue on their github and they said they'll
# fix it soon, but in the meantime I copied the code and made the change myself.
#source("../SimBu/R/dataset.R")
#source("../SimBu/R/simulator.R")

sim_type <- snakemake@wildcards[['pseudobulktype']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Cell types of interest
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells", "ILC",
                "Epithelial cells", "Mast cells", "Plasma cells",
                "Fibroblasts", "Macrophages", "Monocytes", "pDC", "DC")

# Switch simulation names to have sample information and concatenate
reformat_se <- function(overall, specific, i) {
  specific <- as.matrix(assay(specific))
  sample_id <- samples[i]
  sim_count <- dim(specific)[2]
  colnames(specific) <- paste(sample_id, "sim", 1:sim_count, sep = "_")
  if (i == 1) {
    overall <- specific
  } else {
    overall <- cbind(overall, specific)
  }
  
  overall
}

reformat_fractions <- function(overall, specific, i) {
  sample_id <- samples[i]
  sim_count <- nrow(specific)
  rownames(specific) <- paste(sample_id, "sim", 1:sim_count, sep = "_")
  if (i == 1) {
    overall <- specific
  } else {
    overall <- full_join(overall, specific)
    overall[is.na(overall)] <- 0
  }
 
  overall
}

# Create matrix with pseudobulk counts
sim_se <- matrix()

# Create sample by cell type df with true fractions
sim_fractions <- data.frame()

set.seed(1019)
for(i in 1:length(samples)) {
  sample <- samples[i]
  if (sample == "2428"){
    next
  }
  sc_file <- paste(local_data_path, "/sce_objects/", sample, "_labeled.rds", sep = "")
  sce <- readRDS(sc_file)
  sce <- sce[, sce$cellType %in% cell_types]
  colnames(sce) <- sce$Barcode
  rownames(sce) <- rowData(sce)$ID
  
  annotation <- as.data.frame(colData(sce))
  setnames(annotation, "Barcode","ID")
  setnames(annotation, "cellType","cell_type")
  counts <- assay(sce)
  
  data <- dataset(annotation = annotation, count_matrix = counts, name = sample, filter_genes = FALSE)
  if (sim_type == "realistic"){
    sample_sim <- simulate_bulk(data, scenario = "mirror_db",
                                scaling_factor = "NONE",
                                ncells = 2000, 
                                nsamples = 50,
                                balance_even_mirror_scenario = 0.03,
                                remove_bias_in_counts = FALSE,
                                norm_counts = FALSE)
  } else if (sim_type == "even"){
    sample_sim <- simulate_bulk(data, scenario = "even", 
                                scaling_factor = "NONE", 
                                ncells = 2000, 
                                nsamples = 50,
                                remove_bias_in_counts = FALSE,
                                norm_counts = FALSE)
  } else if (sim_type == "weighted"){
    sample_sim <- simulate_bulk(data, scenario = "weighted",
                                scaling_factor = "NONE",
                                ncells = 2000,
                                nsamples = 50,
                                weighted_cell_type = "Epithelial cells",
                                weighted_amount = 0.7,
                                remove_bias_in_counts = FALSE,
                                norm_counts = FALSE)
  } else if (sim_type == "sparse"){
    sample_sim <- simulate_bulk(data, scenario = "mirror_db",
                                scaling_factor = "NONE",
                                ncells = 2000,
                                nsamples = 50, 
                                balance_even_mirror_scenario = 0.03,
                                whitelist = c("Epithelial cells", "T cells", "Fibroblasts",
                                              "Endothelial cells", "Macrophages"),
                                remove_bias_in_counts = FALSE,
                                norm_counts = FALSE)
  }
  sim_se <- reformat_se(sim_se, sample_sim$bulk, i)
  sim_fractions <- reformat_fractions(sim_fractions, sample_sim$cell_fractions, i)
}
rownames(sim_fractions) <- colnames(sim_se)

# Write raw count and cpm normalized data to files
make_files <- function(data, fractions, sim_type) {
  sample_names <- colnames(data)
  cpm_data <- data*1000000/colSums(data)[col(data)]
  rownames(cpm_data) <- rowData(sce)$Symbol
  
  data <- cbind(rownames(data), data)
  colnames(data) <- c("Gene", sample_names)
  
  datafile <- paste(local_data_path, "/deconvolution_input/bulk_data_", sim_type, ".tsv", sep = "")
  write.table(data, datafile, row.names = F, quote = F, sep = "\t")
  
  cpm_data <- cbind(rownames(cpm_data), cpm_data)
  colnames(cpm_data) <- c("Gene", sample_names)
  
  normfile <- paste(local_data_path, "/deconvolution_input/normalized_data_", sim_type, ".tsv", sep = "")
  write.table(cpm_data, normfile, row.names = F, quote = F, sep = "\t")
  
  fractionfile <- paste(local_data_path, "/deconvolution_input/cell_type_fractions_", sim_type, ".tsv", sep = "")
  write.table(fractions, fractionfile, quote = F, sep = "\t")
}

make_files(sim_se, sim_fractions, sim_type)
