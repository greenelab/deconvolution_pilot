# Not only do we have three kinds of bulk RNA-seq, we also have two kinds of
# scRNA-seq: the samples run individually and the pooled samples. Because of
# this, we can use pooled cells as our reference profiles and use the solo
# samples' cell type proportions as a kind of silver standard. We've annotated
# the cell types in the solo samples, and our assumption is that the closer the
# deconvolution results are to the single-cell proportions, the closer the
# deconvolution method gets to the true proportions.

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

true_bulk_types <- c("chunk_ribo", "dissociated_polyA", "dissociated_ribo")

# Cell types of interest
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells",
                "Macrophages", "Mast cells", "Plasma cells", "Epithelial cells",
                "Fibroblasts", "Monocytes", "pDC", "DC", "ILC")

# Load in labeled single cells, melted into a single data frame
melted_sc <- load_melted_sc()

# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% proportions)
melted_results <- subset(melted_results, melted_results$bulk_type %in% true_bulk_types)

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# Filter down to cell types we have single-cell data for
melted_results <- subset(melted_results, melted_results$cell_type %in% cell_types)

# Join single cell proportions to deconvolution results
melted_results <- left_join(melted_results, melted_sc)
melted_results[is.na(melted_results$proportion.sc),]$proportion.sc <- 0

# Plot correlations
plotfile <- paste(plot_path, "/deconvolution_plots/true_bulk_correlations.png", sep = "")
png(plotfile, width = 700)
coors_field <- melted_results %>% group_by(method, bulk_type) %>% summarize(cor = cor(proportion, proportion.sc))
ggplot(coors_field, mapping = aes(x=bulk_type, y=cor, group=method, color=method)) + geom_point() + geom_line()
dev.off()

# Subtract single cell proportions from deconvolution estimates
melted_results$proportion <- melted_results$proportion - melted_results$proportion.sc
melted_results$proportion.sc <- NULL

# Compare accuracy across methods, stratified by cell type and vice versa
plotfile <- paste(plot_path, "/deconvolution_plots/accuracy_by_cell_type.png", sep = "")
png(plotfile, width = 1200)
ggplot(melted_results, aes(x = cell_type, y = proportion, fill = method)) + geom_boxplot() +
	ylab("Estimated proportion - single cell proportion") + xlab("Cell type")
dev.off()

plotfile <- paste(plot_path, "/deconvolution_plots/accuracy_by_method.png", sep = "")
png(plotfile, width = 1200)
ggplot(melted_results, aes(x = method, y = proportion, fill = cell_type)) + geom_boxplot() +
	ylab("Estimated proportion - single cell proportion") + xlab("Cell type")
dev.off()


# Make heatmap of proportion differences for each bulk type separately
chunk_ribo <- subset(melted_results, melted_results$bulk_type == "chunk_ribo")
chunk_ribo$bulk_type <- NULL
chunk_ribo <- as.data.frame(chunk_ribo %>% group_by(cell_type, method) %>% summarize(proportion = mean(proportion)))
chunk_ribo <- reshape(data = chunk_ribo,
                      idvar = "cell_type",
                      timevar = "method",
                      direction = "wide")
rownames(chunk_ribo) <- chunk_ribo$cell_type; chunk_ribo$cell_type <- NULL
setnames(chunk_ribo, c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq"))

plotfile <- paste(plot_path, "/deconvolution_plots/accuracy_heatmap_chunk_ribo.png", sep = "")
png(plotfile, width = 1200)
pheatmap::pheatmap(chunk_ribo, cluster_rows = FALSE, cluster_cols = FALSE,
		   display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks = ppoints(100) * 1.2 - 0.6)
dev.off()

dissociated_ribo <- subset(melted_results, melted_results$bulk_type == "dissociated_ribo")
dissociated_ribo$bulk_type <- NULL
dissociated_ribo <- as.data.frame(dissociated_ribo %>% group_by(cell_type, method) %>% summarize(proportion = mean(proportion)))
dissociated_ribo <- reshape(data = dissociated_ribo,
                      idvar = "cell_type",
                      timevar = "method",
                      direction = "wide")
rownames(dissociated_ribo) <- dissociated_ribo$cell_type; dissociated_ribo$cell_type <- NULL
setnames(dissociated_ribo, c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq"))

plotfile <- paste(plot_path, "/deconvolution_plots/accuracy_heatmap_dissociated_ribo.png", sep = "")
png(plotfile, width = 1200)
pheatmap::pheatmap(dissociated_ribo, cluster_rows = FALSE, cluster_cols = FALSE,
		   display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks = ppoints(100) * 1.2 - 0.6)
dev.off()

dissociated_polyA <- subset(melted_results, melted_results$bulk_type == "dissociated_polyA")
dissociated_polyA$bulk_type <- NULL
dissociated_polyA <- as.data.frame(dissociated_polyA %>% group_by(cell_type, method) %>% summarize(proportion = mean(proportion)))
dissociated_polyA <- reshape(data = dissociated_polyA,
                      idvar = "cell_type",
                      timevar = "method",
                      direction = "wide")
rownames(dissociated_polyA) <- dissociated_polyA$cell_type; dissociated_polyA$cell_type <- NULL
setnames(dissociated_polyA, c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq"))

plotfile <- paste(plot_path, "/deconvolution_plots/accuracy_heatmap_dissociated_polyA.png", sep = "")
png(plotfile, width = 1200)
pheatmap::pheatmap(dissociated_polyA, cluster_rows = FALSE, cluster_cols = FALSE,
		   display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks = ppoints(100) * 1.2 - 0.6)
dev.off()
