# Since we're running three kinds of bulk RNA-seq (one with tumor chunks and
# ribo depletion, one with dissociated cells and ribo depletion, one with
# dissociated cells and poly-A capture), we have three sets of deconvolution
# results for each tumor. We want to see how much proportion estimates vary
# based on the type of bulk data. Our assumption is that if a method gives
# results that are more similar/have a lower variance, they will be more
# robust to technical noise.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
samples <- params$samples

source("evaluation_functions.R")

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq")
scores <- c("consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")

# Cell types of interest
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells",
                "Macrophages", "Mast cells", "Plasma cells", "pDC",
                "Fibroblasts", "Monocytes", "Epithelial cells", "DC")

# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% proportions)

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# For each method x cell type x sample combination, calculate variance across bulk types
results <- melted_results %>% group_by(method, cell_type, variable) %>% summarize(variance = var(proportion))
results <- subset(results, results$cell_type %in% cell_types)

# Compare variance across methods overall
plotfile <- paste(plot_path, "/deconvolution_plots/variance_summary.png", sep = "")
png(plotfile)
ggplot(results, mapping = aes(x = method, y = variance)) + geom_boxplot()
dev.off()

# Compare variance across methods, stratified by cell type and vice versa
plotfile <- paste(plot_path, "/deconvolution_plots/variance_by_method.png", sep = "")
png(plotfile)
ggplot(results, mapping = aes(x = method, y = log(variance), fill = cell_type)) + geom_boxplot()
dev.off()

plotfile <- paste(plot_path, "/deconvolution_plots/variance_by_cell_type.png", sep = "")
png(plotfile)
ggplot(results, mapping = aes(x = cell_type, y = log(variance), fill = method)) + geom_boxplot()
dev.off()
