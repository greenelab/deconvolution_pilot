# Another way of checking if hash demultiplexing vs. genetic demultiplexing
# matters. Since hash demultiplexing with the default cellranger multi
# parameters left out so many cells, we wanted to check if not retrieving
# those cells significantly affects downstream deconvolution results.

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
pseudobulk_types <- c("realistic", "even", "weighted", "sparse")

# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results(demultiplex_default = TRUE)

# Split deconvolution results into demultiplex defaults and original
demultiplexed_types <- grep("demultiplex_default", unique(melted_results$bulk_type), value = T)
demultiplexed <- subset(melted_results, melted_results$bulk_type %in% demultiplexed_types)
melted_results <- subset(melted_results, melted_results$method %in% demultiplexed$method)
original <- subset(melted_results, !melted_results$bulk_type %in% demultiplexed_types)

# Check variance
melted_results$bulk_type <- gsub("_demultiplex_default", "", melted_results$bulk_type)
variance <- melted_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))

plotfile <- paste(plot_path, "/deconvolution_plots/demultiplexing_variance_by_method.png", sep = "")
png(plotfile)
ggplot(variance, mapping = aes(x = method, y = log(variance), fill = cell_type)) + geom_boxplot()
dev.off()

plotfile <- paste(plot_path, "/deconvolution_plots/demultiplexing_variance_by_cell_type.png", sep = "")
png(plotfile)
ggplot(variance, mapping = aes(x = cell_type, y = log(variance), fill = method)) + geom_boxplot()
dev.off()


# Join original and demultiplex results
demultiplexed$bulk_type <- gsub("_demultiplex_default", "", demultiplexed$bulk_type)
setnames(demultiplexed, "proportion", "demultiplexed_proportion")
setnames(original, "proportion", "original_proportion")
deconvolution <- full_join(demultiplexed, original)

# Check correlations
corrs <- deconvolution %>% group_by(method, bulk_type) %>% 
  summarise(cor = cor(demultiplexed_proportion, original_proportion))

plotfile <- paste(plot_path, "/deconvolution_plots/demultiplexing_correlations.png", sep = "")
png(plotfile)
ggplot(corrs, mapping = aes(x=bulk_type, y=cor, group=method, color=method)) + geom_point() +
  geom_line() + xlab("Bulk type") + ylab("Correlation value")
dev.off()

# Check proportion differences
deconvolution$diff <- deconvolution$demultiplexed_proportion - deconvolution$original_proportion

plotfile <- paste(plot_path, "/deconvolution_plots/demultiplexing_accuracy_by_method.png", sep = "")
png(plotfile)
ggplot(deconvolution, mapping = aes(x=method, y=diff, fill = bulk_type)) + geom_boxplot() +
  xlab("Method") + ylab("(Demultiplexed proportion - original proportion)")
dev.off()

plotfile <- paste(plot_path, "/deconvolution_plots/demultiplexing_accuracy_by_bulk_type.png", sep = "")
png(plotfile)
ggplot(deconvolution, mapping = aes(x=bulk_type, y=diff, fill = method)) + geom_boxplot() +
  xlab("Bulk type") + ylab("(Demultiplexed proportion - original proportion)")
dev.off()
