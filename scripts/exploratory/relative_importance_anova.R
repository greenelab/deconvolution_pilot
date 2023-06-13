# A reviewer recommended we run an ANOVA to determine which experimental factors
# make the biggest different on deconvolution results. I had never done an ANOVA
# quite like this one before, particularly given that the various cell type
# estimates are not independent (have a sum-to-one constraint). This script was
# playing around with possible models to workshop what might be appropriate. We
# added a section to the results in the paper discussing the results from the
# chosen model.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(DESeq2)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

source("../evaluation/evaluation_functions.R")

pseudobulk_types <- c("realistic", "even", "weighted", "sparse")
proportions <- c("bayesprism", "bisque", "cibersortx", "epic", "music", "nnls")

# Load real data; we're excluding pseudo-bulk data because cell type proportions
# are very affected by simulation scenario, which isn't relevant to experimental
# design.
melted_results <- load_melted_results(reference_comp = T)
melted_results <- subset(melted_results, !bulk_type %in% pseudobulk_types)
melted_results <- subset(melted_results, method %in% proportions)

# Splitting out bulk type into dissociation status and mRNA enrichment
melted_results$dissociation_status <- ifelse(melted_results$bulk_type=="chunk_ribo", "no", "yes")
melted_results$mRNA_enrichment <- ifelse(melted_results$bulk_type=="dissociated_polyA", "polyA", "rRNA")

# Basic model one
full_model <- aov(proportion ~ method + dissociation_status + mRNA_enrichment + reference, data = melted_results)
summary(full_model)

# Try with just epithelial cells
epithelial <- subset(melted_results, melted_results$cell_type=="Epithelial cells")
epi_model <- aov(proportion ~ method + dissociation_status + mRNA_enrichment + reference, data = epithelial)
summary(epi_model)

# Try with just fibroblasts
fibroblasts <- subset(melted_results, melted_results$cell_type=="Fibroblasts")
fibro_model <- aov(proportion ~ method + dissociation_status + mRNA_enrichment + reference, data = fibroblasts)
summary(fibro_model)

# Try with just endothelial cells
endothelial <- subset(melted_results, melted_results$cell_type=="Endothelial cells")
two_way <- aov(proportion ~ method + dissociation_status + mRNA_enrichment + reference, data = endothelial)
summary(two_way)

# Try with just T cells
tcells <- subset(melted_results, melted_results$cell_type=="T cells")
two_way <- aov(proportion ~ method + dissociation_status + mRNA_enrichment + reference, data = tcells)
summary(two_way)

# Try with just BayesPrism
bp <- subset(melted_results, melted_results$method == "bayesprism")
bayes <- aov(proportion ~ dissociation_status + mRNA_enrichment + reference, data = bp)
summary(bayes)

# Basic model with sample included
sample_model <- aov(proportion ~ variable + method + dissociation_status + mRNA_enrichment, data = epithelial)
summary(sample_model)