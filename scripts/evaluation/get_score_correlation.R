# Some deconvolution methods don't return proportions, but instead cell type
# scores that are meant to be compared across samples and don't have a real
# unit. Because of this, sometimes people doing benchmarking for these methods
# will compare the correlations between the scores and known values. For this,
# we'll look at the correlation between scores and single-cell proportions.

suppressPackageStartupMessages({
	library(data.table)
	library(SingleCellExperiment)
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
cell_types <- c("NK cells", "CD4 T cells", "CD8 T cells", "Regulatory T cells",
		"B cells", "Endothelial cells", "Macrophages", "Mast cells",
		"Plasma cells", "Fibroblasts", "Monocytes", "pDC", "DC")

# Load in labeled single cells, melted into a single data frame
melted_sc <- load_melted_sc(granular = TRUE)

# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% scores)
setnames(melted_results, "variable", "sample")

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# Plot correlations
for (i in 1:length(scores)){
	methodname <- scores[i]
	x <- subset(melted_results, melted_results$method == methodname)
	x <- subset(x, x$sample != "2428")
	mutual_cell_types <- intersect(x$cell_type, melted_sc$cell_type)
	x <- full_join(x, melted_sc)
	x <- subset(x, x$cell_type %in% mutual_cell_types)

	# Mark empty categories in single cell as proportions of 0
	x[is.na(x$proportion.sc), ]$proportion.sc <- 0

	plotfile <- paste(plot_path, "/deconvolution_plots/", method_name, "_correlation.png", sep = "")
	png(plotfile)
	title <- paste(methodname, ", r = ", cor(x$score, x$proportion.sc), ", rho = ",
		       cor(x$score, x$proportion.sc, method = "spearman"), sep = "")
	ggplot(x, mapping = aes(x = score, y = proportion.sc, color = bulk_type)) +
		geom_point() + ggtitle(title)
	dev.off()
}
