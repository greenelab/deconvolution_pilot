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
figure_path <- params$figure_path
samples <- params$samples

source("figure_utils.R")
source("../evaluation/evaluation_functions.R")

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("bayesprism", "bisque", "cibersortx", "epic", "music", "nnls")
scores <- c("consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")
single_cell_methods <- c("bayesprism","cibersortx","bisque","music","nnls")

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

pA <- ggplot(results) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across bulk types)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes)


# Load all deconvolution results across reference profiles
melted_results <- load_melted_results(reference_comp = TRUE)
reference_methods <- unique(subset(melted_results, melted_results$reference=="hashing")$method)
melted_results <- subset(melted_results, melted_results$method %in% reference_methods)

# Split into only non-simulated reference profiles and plot variance
nonsim_results <- subset(melted_results, melted_results$reference %in%
                           c("hashing", "genetic"))
results <- nonsim_results %>% group_by(method, cell_type, variable, bulk_type) %>% summarize(variance = var(proportion))

pB <- ggplot(results) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across reference profiles)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes)

# Join hashing and genetic reference results
genetic_results <- subset(nonsim_results, reference=="genetic")
genetic_results$reference <- NULL; setnames(genetic_results, "proportion", "genetic_proportion")
hashing_results <- subset(nonsim_results, reference=="hashing")
hashing_results$reference <- NULL; setnames(hashing_results, "proportion", "hashing_proportion")
deconvolution <- full_join(genetic_results, hashing_results)

# Get RMSE
deconvolution$prop_diff <- deconvolution$genetic_proportion - deconvolution$hashing_proportion
deconvolution$prop_diff_sq <- deconvolution$prop_diff ^ 2
sum_sqs <- deconvolution %>% group_by(method, bulk_type) %>% summarize(rmse = sqrt(mean(prop_diff_sq)))

# Check RMSE across bulk data types
sum_sqs$bulk_type <- recode(sum_sqs$bulk_type,
		       "chunk_ribo" = "rRNA- Chunk",
		       "dissociated_ribo" = "rRNA- Dissociated",
		       "dissociated_polyA" = "polyA+ Dissociated")
sum_sqs$bulk_type <- factor(sum_sqs$bulk_type, levels = c("rRNA- Chunk", "rRNA- Dissociated",
							 "polyA+ Dissociated", "even",
							 "realistic", "sparse", "weighted"))

pC <- ggplot(sum_sqs, mapping = aes(x = bulk_type, y = rmse, group = method, color = method)) + geom_point() +
  geom_line() + xlab("Bulk type") + ylab("RMSE across reference profiles") +
  scale_color_manual(name = "Method", values = colors_methods,
                     limits = c("bayesprism", "bisque","cibersortx","music","nnls"))


# Check proportion differences
deconvolution$bulk_type <- recode(deconvolution$bulk_type,
			       "chunk_ribo" = "rRNA- Chunk",
			       "dissociated_ribo" = "rRNA- Dissociated",
			       "dissociated_polyA" = "polyA+ Dissociated")
deconvolution$bulk_type <- factor(deconvolution$bulk_type, levels = c("rRNA- Chunk", "rRNA- Dissociated",
							 "polyA+ Dissociated", "even",
							 "realistic", "sparse", "weighted"))


pD <- ggplot(deconvolution) + geom_boxplot(mapping = aes(x = method, y = prop_diff, fill = bulk_type, color = bulk_type)) +
  geom_boxplot(mapping = aes(x = method, y = prop_diff, fill = bulk_type), outlier.colour = NA) +
  xlab("Method") + ylab("Proportion (genetic reference - hashing reference)") +
  scale_color_manual(name = "Bulk type", values = colors_bulktypes) +
  scale_fill_manual(name = "Bulk type", values = colors_bulktypes)


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
robustness$average_var <- log(robustness$average_var)

# Load pseudobulk true fractions
melted_fractions <- load_pseudobulk_fractions(pseudobulk_types)
setnames(melted_fractions, "proportion", "true_proportion")

# Filter down to pseudobulk data
melted_pseudo_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
melted_pseudo_results <- subset(melted_pseudo_results, melted_pseudo_results$cell_type %in% melted_pseudo_results$cell_type)

# Join true proportions to deconvolution results
melted_pseudo_results <- left_join(melted_pseudo_results, melted_fractions)
melted_pseudo_results[is.na(melted_pseudo_results$true_proportion), ]$true_proportion <- 0

# Get RMSE
melted_pseudo_results$prop_diff <- melted_pseudo_results$proportion - melted_pseudo_results$true_proportion
melted_pseudo_results$prop_diff_sq <- melted_pseudo_results$prop_diff ^ 2
sum_sqs <- melted_pseudo_results %>% group_by(method) %>% summarize(rmse = sqrt(mean(prop_diff_sq)))

# Load real data proportions
melted_sc <- load_melted_sc()
melted_real_results <- left_join(melted_real_results, melted_sc)
melted_real_results[is.na(melted_real_results$proportion.sc), ]$proportion.sc <- 0

# Get average sum of least squares
melted_real_results$prop_diff <- melted_real_results$proportion - melted_real_results$proportion.sc
melted_real_results$prop_diff_sq <- melted_real_results$prop_diff ^ 2
real_sum_sqs <- melted_real_results %>% group_by(method) %>% summarize(real_rmse = sqrt(mean(prop_diff_sq)))

# Plot accuracy vs robustness
total <- full_join(robustness, sum_sqs)
pE <- ggplot(total, mapping = aes(x = average_var, y = rmse, color = method)) +
  geom_point(aes(size = 10)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (log(variance) across bulk types)") +
  ylab("Accuracy (RMSE with pseudo-bulk proportions)") +
  guides(color = "none", size = "none") +
  scale_color_manual(name = "Method", values = colors_methods)

total <- full_join(robustness, real_sum_sqs)
pF <- ggplot(total, mapping = aes(x = average_var, y = real_rmse, color = method)) +
  geom_point(aes(size = 10)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  geom_label_repel(aes(label = method, size = NULL)) +
  xlab("Robustness (log(variance) across bulk types)") +
  ylab("Accuracy (RMSE with single cell proportions)") +
  guides(color = "none", size = "none") +
  scale_color_manual(name = "Method", values = colors_methods)

pdf(paste(figure_path, "figure6.pdf", sep = "/"), width = 16, height = 16, family = "sans")
top <- pA + pB
middle <- (pC + theme(axis.title.x = element_text(margin = margin(t = -30, unit = "pt")))) + (pD + theme(axis.title.x = element_text(margin = margin(t = -30, unit = "pt")))) + plot_layout(ncol = 2, width = c(2, 4))
bottom <- pE + plot_spacer() + pF + plot_layout(ncol = 3, widths = c(3, 1, 3))
top / middle / bottom + plot_annotation(tag_levels = "A")
dev.off()




## Variance of deconvolution results with smaller and smaller reference profiles
melted_results <- load_melted_results(reference_comp = TRUE)
melted_results <- rename_cell_types(melted_results)
melted_results <- subset(melted_results, melted_results$method %in% single_cell_methods)

## Split up results into bulk and pseudo and identify all possible combos
melted_real_results <- subset(melted_results, !melted_results$bulk_type %in% pseudobulk_types)
all_real_combos <- melted_real_results %>%
  tidyr::expand(cell_type, variable, bulk_type, reference, method) 
melted_real_results <- melted_real_results %>% right_join(all_real_combos)

melted_pseudo_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
all_pseudo_combos <- melted_pseudo_results %>%
  tidyr::expand(cell_type, variable, bulk_type, reference, method) 
melted_pseudo_results <- melted_pseudo_results %>% right_join(all_pseudo_combos)

melted_results <- rbind(melted_real_results, melted_pseudo_results)
melted_results <- melted_results %>% mutate(proportion = ifelse(is.na(proportion), 0, proportion))

variance3 <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000")) %>%
  group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
sA <- ggplot(variance3) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across reference profiles)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes) + ylim(-45, -3) +
  ggtitle("Genetic, hashing, Sim2000")
variance3 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

variance4 <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000","sim1000")) %>% 
  group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
sB <- ggplot(variance4) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across reference profiles)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes) + ylim(-45, -3) +
  ggtitle("Adding Sim1000")
variance4 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

variance5 <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000","sim1000","sim500")) %>%
  group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
sC <- ggplot(variance5) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across reference profiles)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes) + ylim(-45, -3) +
  ggtitle("Adding Sim500")
variance5 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

variance6 <- melted_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
sD <- ggplot(variance6) + geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type, color = cell_type)) +
  geom_boxplot(mapping = aes(x = method, y = log(variance), fill = cell_type), outlier.colour = NA) +
  xlab("Method") + ylab("log(variance across reference profiles)") +
  scale_color_manual(name = "Cell type", values = colors_celltypes) +
  scale_fill_manual(name = "Cell type", values = colors_celltypes) + ylim(-45, -3) +
  ggtitle("Adding Sim200")
variance6 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

pdf(paste(figure_path, "suppfig6.pdf", sep = "/"), width = 16, height = 12, family = "sans")
sA + sB + sC + sD + plot_annotation(tag_levels = "A")
dev.off()
