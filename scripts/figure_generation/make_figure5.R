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
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

source("figure_utils.R")
source("../evaluation/evaluation_functions.R")

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("bayesprism", "bisque", "cibersortx", "epic", "music", "nnls")
scores <- c("consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")

# Plot true proportions for pseudobulk samples
pseudobulk_types <- c("realistic", "even", "weighted", "sparse")
melted_fractions <- load_pseudobulk_fractions(pseudobulk_types)

# Filter down to one sample for easier visualization
demo_fractions <- melted_fractions[grep("2380", melted_fractions$sample),]

pA <-ggplot(demo_fractions, mapping = aes(x=sample, y=proportion, fill=cell_type)) +
  geom_bar(stat = "identity") +
  #theme(legend.position = "bottom") + 
  facet_wrap("~bulk_type", ncol = 4) +
  scale_fill_manual(values = colors_celltypes) + 
  labs(fill = "Cell type") +
  xlab("Simulated sample") + ylab("Proportion") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) #remove x axis ticks



# Load all pseudobulk deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% proportions)
setnames(melted_results, "variable", "sample")

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

# Split into real bulk and pseudobulk
melted_pseudo_results <- subset(melted_results, melted_results$bulk_type %in% pseudobulk_types)
melted_real_results <- subset(melted_results, !melted_results$bulk_type %in% pseudobulk_types)

# Filter down to cell types we have single-cell data for
melted_pseudo_results <- subset(melted_pseudo_results, melted_pseudo_results$cell_type %in% melted_fractions$cell_type)

# Join true proportions to deconvolution results
setnames(melted_fractions, "proportion", "true_proportion")
melted_pseudo_results <- left_join(melted_pseudo_results, melted_fractions)
melted_pseudo_results[is.na(melted_pseudo_results$true_proportion),]$true_proportion <- 0

# Plot correlations
correlations <- melted_pseudo_results %>%
  group_by(method, bulk_type) %>%
  summarize(cor = cor(proportion, true_proportion))
pB <- ggplot(correlations, mapping = aes(x=bulk_type, y=cor, group=method, color=method)) +
  geom_point() + geom_line() +
  scale_color_manual(name = "Method", values = colors_methods) +
  xlab("Pseudo-bulk simulation type") +
  ylab("Correlation with true pseudo-bulk proportions") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

# Subtract true proportions from deconvolution estimates
melted_pseudo_results$proportion <- melted_pseudo_results$proportion - melted_pseudo_results$true_proportion
melted_pseudo_results$true_proportion <- NULL

# Plot accuracy by cell type
pC <- ggplot(melted_pseudo_results) + geom_boxplot(mapping = aes(x = cell_type, y= proportion, fill = method, color = method)) +
  geom_boxplot(mapping = aes(x = cell_type, y = proportion, fill = method), outlier.colour = NA) +
  scale_color_manual(name = "Method", values = colors_methods) +
  scale_fill_manual(name = "Method", values = colors_methods) +
  ylab("Estimated proportion - true proportion") + xlab("Cell type")




# Load in labeled single cells, melted into a single data frame
melted_sc <- load_melted_sc()

# Filter down to cell types we have single-cell data for
melted_real_results <- subset(melted_real_results, melted_real_results$cell_type %in% melted_sc$cell_type)

# Join single cell proportions to deconvolution results
melted_real_results <- left_join(melted_real_results, melted_sc)
melted_real_results[is.na(melted_real_results$proportion.sc),]$proportion.sc <- 0

# Subtract single cell proportions from deconvolution estimates
melted_real_results$proportion <- melted_real_results$proportion - melted_real_results$proportion.sc
melted_real_results$proportion.sc <- NULL

# Plot true accuracy by cell type
pD <- ggplot(melted_real_results) + geom_boxplot(mapping = aes(x = cell_type, y= proportion, fill = method, color = method)) +
  geom_boxplot(mapping = aes(x = cell_type, y = proportion, fill = method), outlier.colour = NA) +
  scale_color_manual(name = "Method", values = colors_methods) +
  scale_fill_manual(name = "Method", values = colors_methods) +
  ylab("Estimated proportion - single cell proportion") + xlab("Cell type")



# Load MCPcounter genes
mcpcounter <- fread(paste(local_data_path, "miscellaneous", "MCPcounter.txt", sep = "/"))

# Load real bulk vs. pseudobulk DESeq2 object
deseq_path <- paste(local_data_path, "deseq2_output", sep = "/")
dds <- readRDS(paste(deseq_path, "polyA_vs_pseudo_data.rds", sep = "/"))

res <- as.data.frame(results(dds))
res$gene <- rownames(res)
setnames(mcpcounter, "HUGO symbols", "gene")
setnames(mcpcounter, "Cell population", "cell_type")

mcpcounter <- inner_join(mcpcounter, res)

# Make volcano plot
pE <- ggplot(mcpcounter, mapping = aes(x = log2FoldChange, y = -log10(padj), color = cell_type)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  xlab("log2 fold change") + ylab("-log10 adjusted p-value") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  labs(color = "Cell type")

pdf("../../figures/figure5.pdf", width = 20, height = 12, family = "sans")
top <- pA + pB + plot_layout(ncol = 2, width = c(5, 2))
bottom <- pC + pD + pE + plot_layout(ncol = 3, width = c(3, 3, 2))
top / bottom + plot_annotation(tag_levels = "A")
dev.off()


make_proportion_heatmap <- function(melted_results, bt) {
  melted_results <- subset(melted_results, melted_results$bulk_type == bt)
  melted_results$bulk_type <- NULL
  melted_results <- as.data.frame(melted_results %>% group_by(cell_type, method) %>% summarize(proportion = mean(proportion)))
  melted_results$proportion <- round(melted_results$proportion, digits = 2)

  melted_results$cell_type <- factor(melted_results$cell_type, 
                                     levels = sort(unique(melted_results$cell_type),
                                                   decreasing = TRUE))
  
  # Reshape and then melt to replace NAs with 0s
  reshaped_heatmap <- reshape(data=melted_results, 
                       idvar="cell_type",
                       v.names="proportion",
                       timevar= "method",
                       direction="wide")
  colnames(reshaped_heatmap) <- gsub("proportion.", "", colnames(reshaped_heatmap))
  
  remelt <- melt(reshaped_heatmap)
  setnames(remelt, c("cell_type", "method", "proportion"))
  
  remelt
}



# Make heatmap of proportion differences for each bulk type separately
chunk_ribo <- make_proportion_heatmap(melted_real_results, "chunk_ribo")
sA <- ggplot(chunk_ribo, aes(x=method, y=cell_type, fill=proportion)) +
  geom_raster() + geom_text(aes(label = proportion), size = 6) +
scale_fill_gradientn(colors = heatmap_scale_2d, limits = c(-0.7,0.7), na.value = "#DDDDDD") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "None")

dissociated_ribo <- make_proportion_heatmap(melted_real_results, "dissociated_ribo")
sB <- ggplot(dissociated_ribo, aes(x=method, y=cell_type, fill=proportion)) +
  geom_raster() + geom_text(aes(label = proportion), size = 6) +
  scale_fill_gradientn(colors = heatmap_scale_2d, limits = c(-0.7,0.7), na.value = "#DDDDDD") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "None")

dissociated_polyA <- make_proportion_heatmap(melted_real_results, "dissociated_polyA")
sC <- ggplot(dissociated_polyA, aes(x=method, y=cell_type, fill=proportion)) +
  geom_raster() + geom_text(aes(label = proportion), size = 6) +
  scale_fill_gradientn(colors = heatmap_scale_2d, limits = c(-0.7,0.7), na.value = "#DDDDDD") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "None")

pdf("../../figures/suppfig5.pdf", width = 24, height = 12, family = "sans")
sA + sB + sC +
  plot_annotation(tag_levels = "A")
dev.off()
