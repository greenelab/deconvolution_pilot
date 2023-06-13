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
proportions <- c("bayesprism", "bisque", "cibersortx", "epic", "music", "nnls")
scores <- c("quantiseq", "abis", "consensus_tme", "immucellai", "mcpcounter", "timer", "xcell")

# Cell types of interest
cell_types <- c("NK cells", "CD4 T cells", "CD8 T cells", "Regulatory T cells",
                "B cells", "Endothelial cells", "Macrophages", "Mast cells",
                "Plasma cells", "Fibroblasts", "Monocytes", "pDC", "DC")
pseudobulk_types <- c("realistic", "even", "weighted", "sparse")


# Load in labeled single cells, melted into a single data frame
melted_sc <- load_melted_sc(granular = TRUE)


# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results()
melted_results <- subset(melted_results, melted_results$method %in% scores &
                             !melted_results$bulk_type %in% pseudobulk_types) 

# Unify cell type nomenclature across methods
melted_results <- rename_cell_types(melted_results)

setnames(melted_results,"variable", "sample")

# Plot correlations
plot_score <- function(methodname){
    x <- subset(melted_results, melted_results$method == methodname)
    x <- subset(x, x$sample != "2428")
    mutual_cell_types <- intersect(x$cell_type, melted_sc$cell_type)
    x <- full_join(x, melted_sc)
    x <- subset(x, x$cell_type %in% mutual_cell_types)
    
    # Mark empty categories in single cell as proportions of 0
    x[is.na(x$proportion.sc), ]$proportion.sc <- 0
    
    title <- paste(methodname, ", r = ", cor(x$proportion, x$proportion.sc), ", rho = ",
                   cor(x$proportion, x$proportion.sc, method = "spearman"), sep = "")
    qTemp <- ggplot(x, mapping = aes(x = proportion, y = proportion.sc, color = bulk_type)) +
        geom_point() + ggtitle(title)
    
    qTemp
}

qA <- plot_score("abis")
qB <- plot_score("consensus_tme")
qC <- plot_score("immucellai")
qD <- plot_score("mcpcounter")
qE <- plot_score("quantiseq")
qF <- plot_score("timer")
qG <- plot_score("xcell")


pdf("../../figures/suppfig7.pdf", width = 12, height = 6, family = "sans")
qA + qB + qC + qD + qE + qF + qG + 
    plot_layout(nrow = 2) + 
    plot_annotation(tag_levels = "A")
dev.off()
