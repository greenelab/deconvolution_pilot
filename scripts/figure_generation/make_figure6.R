suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(scater)
  library(miQC)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(ggpubr)
  library(yaml)
  library(DESeq2)
  library(rtracklayer)
  library(ggrepel)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

source("figure_utils.R")
source("../evaluation/evaluation_functions.R")

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("abis", "bayesprism", "bisque", "cibersortx", "epic", "music", "nnls", "quantiseq")
scores <- c("consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")

# Cell types of interest
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells",
                "Macrophages", "Mast cells", "Plasma cells", "pDC", "ILC",
                "Fibroblasts", "Monocytes", "Epithelial cells", "DC")
pseudobulk_types <- c("realistic", "even", "weighted", "sparse")


# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% proportions &
                           !melted_results$bulk_type %in% pseudobulk_types) 

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# For each method x cell type x sample combination, calculate variance across bulk types
results <- melted_results %>% group_by(method, cell_type, variable) %>% summarize(variance = var(proportion))
results <- subset(results, results$cell_type %in% cell_types)

pA <- ggplot(results, mapping = aes(x = method, y = log(variance), color = cell_type)) + geom_boxplot()



# Split deconvolution results into demultiplex defaults and original
melted_results <- load_melted_results(demultiplex_default = TRUE)
demultiplexed_types <- grep("demultiplex_default", unique(melted_results$bulk_type), value = T)
demultiplexed <- subset(melted_results, melted_results$bulk_type %in% demultiplexed_types)
melted_results <- subset(melted_results, melted_results$method %in% demultiplexed$method)
original <- subset(melted_results, !melted_results$bulk_type %in% demultiplexed_types)

# Check variance
melted_results$bulk_type <- gsub("_demultiplex_default", "", melted_results$bulk_type)
variance <- melted_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))

pB <- ggplot(variance, mapping = aes(x = method, y = log(variance), color = cell_type)) + geom_boxplot()

# Join original and demultiplex results
demultiplexed$bulk_type <- gsub("_demultiplex_default", "", demultiplexed$bulk_type)
setnames(demultiplexed, "proportion", "demultiplexed_proportion")
setnames(original, "proportion", "original_proportion")
deconvolution <- full_join(demultiplexed, original)

# Check correlations
corrs <- deconvolution %>% group_by(method, bulk_type) %>% 
  summarize(cor = cor(demultiplexed_proportion, original_proportion))
pC <- ggplot(corrs, mapping = aes(x=bulk_type, y=cor, group=method, color=method)) + geom_point() +
  geom_line() + xlab("Bulk type") + ylab("Correlation value")

# Check proportion differences
deconvolution$diff <- deconvolution$demultiplexed_proportion - deconvolution$original_proportion
pD <- ggplot(deconvolution, mapping = aes(x=method, y=diff, color = bulk_type)) + geom_boxplot() +
  xlab("Method") + ylab("(Demultiplexed proportion - original proportion)")




# Load all pseudobulk deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% proportions)
setnames(melted_results, "variable", "sample")

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)
melted_results <- subset(melted_results, melted_results$cell_type %in% cell_types)

# Split into real bulk and pseudobulk
melted_pseudo_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
melted_real_results <- subset(melted_results, !melted_results$bulk_type %in% pseudobulk_types)

# For each method x cell type x sample combination, calculate variance across bulk types
results <- melted_real_results %>% group_by(method, cell_type, sample) %>% summarize(variance = var(proportion))
results <- subset(results, results$cell_type %in% cell_types)

robustness <- results %>% group_by(method) %>% summarize(average_var = mean(variance))
robustness$average_var <- -log(robustness$average_var)

# Load pseudobulk true fractions
melted_fractions <- load_pseudobulk_fractions(pseudobulk_types)
setnames(melted_fractions, "proportion", "true_proportion")

# Filter down to pseudobulk data
melted_pseudo_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
melted_pseudo_results <- subset(melted_pseudo_results, melted_pseudo_results$cell_type %in% melted_pseudo_results$cell_type)

# Join true proportions to deconvolution results
melted_pseudo_results <- left_join(melted_pseudo_results, melted_fractions)
melted_pseudo_results[is.na(melted_pseudo_results$true_proportion), ]$true_proportion <- 0

# Get correlations for pseudobulk data
coors <- melted_pseudo_results %>% group_by(method) %>% summarize(cor = cor(proportion, true_proportion))

# Get average sum of least squares
melted_pseudo_results$prop_diff <- melted_pseudo_results$proportion - melted_pseudo_results$true_proportion
melted_pseudo_results$prop_diff_sq <- melted_pseudo_results$prop_diff ^ 2
sum_sqs <- melted_pseudo_results %>% group_by(method) %>% summarize(rmse = sqrt(mean(prop_diff_sq)))

# Load real data proportions
melted_sc <- load_melted_sc()
melted_real_results <- left_join(melted_real_results, melted_sc)
melted_real_results[is.na(melted_real_results$proportion.sc), ]$proportion.sc <- 0

# Get correlations for real data
real_coors <- melted_real_results %>% group_by(method) %>% summarize(real_cor = cor(proportion, proportion.sc))

# Get average sum of least squares
melted_real_results$prop_diff <- melted_real_results$proportion - melted_real_results$proportion.sc
melted_real_results$prop_diff_sq <- melted_real_results$prop_diff ^ 2
real_sum_sqs <- melted_real_results %>% group_by(method) %>% summarize(real_rmse = sqrt(mean(prop_diff_sq)))

## Completeness (if the method gives results for all cell types of interest)
completeness_results <- subset(melted_real_results, melted_real_results$sample != "2428")
completeness <- completeness_results %>% group_by(method, bulk_type, sample)%>% summarize(possible_proportions = sum(proportion.sc)) %>%
  group_by(method) %>% summarize(pct_explained = mean(possible_proportions))
completeness$immune_only <- completeness$method %in% c("quantiseq", "abis")


# Plot accuracy vs robustness
total <- full_join(robustness, coors)
total <- full_join(total, completeness)
pE <- ggplot(total, mapping = aes(x = average_var, y = cor, color = method, shape = immune_only)) +
  geom_point(aes(size = 10)) + theme(text = element_text(size = 14)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (-log(variance) across protocols)") +
  ylab("Accuracy (correlation with pseudobulk fractions)") +
  guides(color = "none", size = "none")

total <- full_join(robustness, real_coors)
total <- full_join(total, completeness)
pF <- ggplot(total, mapping = aes(x = average_var, y = real_cor, color = method, shape = immune_only)) +
  geom_point(aes(size = 10)) + theme(text = element_text(size = 14)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (-log(variance) across protocols)") +
  ylab("Accuracy (correlation with single cell fractions)") +
  guides(color = "none", size = "none")




pdf("../../figures/figure6.pdf", width = 12, height = 14, family = "sans")
pA + pB + pC + pD + pE + pF +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")
dev.off()





## Accuracy vs. robustness plots using RMSE

total <- full_join(robustness, sum_sqs)
total <- full_join(total, completeness)
qA <- ggplot(total, mapping = aes(x = average_var, y = rmse, color = method, shape = immune_only)) +
  geom_point(aes(size = 10)) + theme(text = element_text(size = 14)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (-log(variance) across protocols)") +
  ylab("Accuracy (RMSE compared to pseudobulk fractions)") +
  guides(color = "none", size = "none")

total <- full_join(robustness, real_sum_sqs)
total <- full_join(total, completeness)
qB <- ggplot(total, mapping = aes(x = average_var, y = real_rmse, color = method, shape = immune_only)) +
  geom_point(aes(size = 10)) + theme(text = element_text(size = 14)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (-log(variance) across protocols)") +
  ylab("Accuracy (RMSE compared to single cell fractions)") +
  guides(color = "none", size = "none")

pdf("../../figures/suppfig6.pdf", width = 12, height = 6, family = "sans")
qA + qB + plot_annotation(tag_levels = "A")
dev.off()
