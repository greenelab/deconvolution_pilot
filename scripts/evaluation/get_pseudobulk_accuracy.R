# With pseudobulk data, it's possible to know the "true" cell type proportions
# and make direct comparisons. Here we will measure the relationships between
# deconvolution results and true proportions in two ways: by computing their
# correlations (what's normally done for benchmarks) and by plotting the 
# deconvolution results minus the true proportions by cell type. 

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

source("evaluation_functions.R")

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq")
scores <- c("consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")

pseudobulk_types <- c("realistic", "even", "weighted", "sparse")

# Load all pseudobulk deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
melted_results <- subset(melted_results, melted_results$method %in% proportions)
setnames(melted_results, "variable", "sample")

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# Load pseudobulk true fractions
melted_fractions <- load_pseudobulk_fractions(pseudobulk_types)
setnames(melted_fractions, "proportion", "true_proportion")

# Filter down to cell types we have single-cell data for
melted_results <- subset(melted_results, melted_results$cell_type %in% melted_fractions$cell_type)

# Join true proportions to deconvolution results
melted_results <- left_join(melted_results, melted_fractions)
melted_results[is.na(melted_results$true_proportion),]$true_proportion <- 0

# Plot correlations
plotfile <- paste(plot_path, "/deconvolution_plots/pseudobulk_correlations.png", sep = "")
png(plotfile, width = 700)
coors_field <- melted_results %>% group_by(method, bulk_type) %>% summarise(cor = cor(proportion, true_proportion))
ggplot(coors_field, mapping = aes(x=bulk_type, y=cor, group=method, color=method)) + geom_point() + geom_line()
dev.off()

# Subtract true proportions from deconvolution estimates
melted_results$proportion <- melted_results$proportion - melted_results$true_proportion
melted_results$true_proportion <- NULL

# Compare accuracy across methods, stratified by cell type and vice versa
plotfile <- paste(plot_path, "/deconvolution_plots/pseudobulk_accuracy_by_cell_type.png", sep = "")
png(plotfile, width = 700)
ggplot(melted_results, aes(x = cell_type, y = proportion, fill = method)) + geom_boxplot() +
  ylab("Estimated proportion - true proportion") + xlab("Cell type")
dev.off()

plotfile <- paste(plot_path, "/deconvolution_plots/pseudobulk_accuracy_by_method.png", sep = "")
png(plotfile, width = 700)
ggplot(melted_results, aes(x = method, y = proportion, fill = cell_type)) + geom_boxplot() +
  ylab("Estimated proportion - true proportion") + xlab("Cell type")
dev.off()

# Split out by type of pseudobulk simulation
for(i in 1:length(pseudobulk_types)){
  scenario <- pseudobulk_types[i]
  scenario_data <- subset(melted_results, melted_results$bulk_type==scenario)
  
  plotfile <- paste(plot_path, "/deconvolution_plots/pseudobulk_", scenario, "_accuracy_by_cell_type.png", sep = "")
  png(plotfile, width = 700)
  ggplot(scenario_data, aes(x = cell_type, y = proportion, fill = method)) + geom_boxplot() +
    ylab("Estimated proportion - true proportion") + xlab("Cell type") + ggtitle(scenario)
  dev.off()
  
  plotfile <- paste(plot_path, "/deconvolution_plots/pseudobulk_", scenario, "_accuracy_by_method.png", sep = "")
  png(plotfile, width = 700)
  ggplot(scenario_data, aes(x = method, y = proportion, fill = cell_type)) + geom_boxplot() +
    ylab("Estimated proportion - true proportion") + xlab("Method") + ggtitle(scenario)
  dev.off()
}

