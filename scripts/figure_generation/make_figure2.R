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
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

source("figure_utils.R")

load_and_filter_datasets <- function(sample_id) {
  # Load the datasets
  hashing <- fread(paste(data_path, "pooled_tumors", sample_id,
                             "Cellranger/outs/multi/multiplexing_analysis",
                             "assignment_confidence_table.csv", sep = "/"))
  genetic <- fread(paste(data_path, "pooled_tumors", sample_id, 
                         "vireo/chunk_ribo/donor_ids.tsv", sep = "/"))
  sce <- readRDS(paste(local_data_path,"/sce_objects/",
                       sample_id, ".rds", sep = ""))
  
  # Cut off "anti-human" part of hashtag labels for legend labels
  hashing$Assignment <- gsub("anti-human_","",hashing$Assignment)
  
  # Rename assignment column to match future label convention
  setnames(hashing, "Assignment", "assignment_0.90")
  
  # Rename barcode columns to match
  setnames(genetic, "cell", "Barcodes")
  
  # Set assignment columns to factors, so they'll order properly in confusion matrix
  hashing$assignment_0.90 <- as.factor(hashing$assignment_0.90)
  genetic$donor_id <- as.factor(genetic$donor_id)
  
  # Filter hashing and vireo to cells in sce object, all others failed miQC filtering
  hashing <- subset(hashing, hashing$Barcodes %in% sce$Barcode)
  genetic <- subset(genetic, genetic$Barcodes %in% sce$Barcode)
  
  list(hashing, genetic, sce)
}

load_cell_types <- function(sample_id, hashing) {
  # Load celltypist assignments
  ct <- fread(paste(local_data_path, "/celltypist_output/",
                    sample_id, "_predicted_labels.csv", sep = ""))
  
  # Consolidate based on hierarchical structure from CellTypist encyclopedia
  label_table <- fread(paste(local_data_path,
                             "/celltypist_output/simplified_labels.tsv",
                             sep = ""))
  setnames(label_table, "Original", "majority_voting")
  ct <- left_join(ct, label_table)
  
  # Combine cell type and sample assignments
  setnames(ct, "V1", "Barcodes")
  ct <- left_join(ct, hashing)
  
  ct
}

set_new_threshold <- function(assignment_table,
                              hashtag_prob, 
                              other_hashtag_prob = 0.001) {
  assignment_table[assignment_table$assignment_0.90 == "Unassigned" &
                     assignment_table$`anti-human_Hashtag1` >= hashtag_prob & 
                     assignment_table$`anti-human_Hashtag2` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag3` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag4` <= other_hashtag_prob, ]$assignment_0.90 = "Hashtag1"
  assignment_table[assignment_table$assignment_0.90 == "Unassigned" &
                     assignment_table$`anti-human_Hashtag1` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag2` >= hashtag_prob &
                     assignment_table$`anti-human_Hashtag3` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag4` <= other_hashtag_prob,
                   ]$assignment_0.90 = "Hashtag2"
  assignment_table[assignment_table$assignment_0.90 == "Unassigned" &
                     assignment_table$`anti-human_Hashtag1` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag2` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag3` >= hashtag_prob &
                     assignment_table$`anti-human_Hashtag4` <= other_hashtag_prob,
                   ]$assignment_0.90 = "Hashtag3"
  assignment_table[assignment_table$assignment_0.90 == "Unassigned" &
                     assignment_table$`anti-human_Hashtag1` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag2` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag3` <= other_hashtag_prob &
                     assignment_table$`anti-human_Hashtag4` >= hashtag_prob,
                   ]$assignment_0.90 = "Hashtag4"
  
  assignment_table$assignment_0.90
}

make_overlap_matrix <- function(hashing,
                                genetic,
                                hashtag_prob) {
  master_list <- full_join(hashing, genetic)

  hash_assignment_label <- paste("assignment", hashtag_prob, sep = "_")
  groups <- c(hash_assignment_label, "donor_id")
  master_list <- subset(master_list, select = groups)
  grouped_data <- master_list %>% 
    group_by(across(all_of(groups))) %>%
    summarize(count = n()) %>% as.data.frame()
  
  # Reshape and then melt to replace NAs with 0s
  confusion <- reshape(data=grouped_data, 
                       idvar="donor_id",
                       v.names="count",
                       timevar=hash_assignment_label,
                       direction="wide")
  colnames(confusion) <- gsub("count.", "", colnames(confusion))
  confusion[is.na(confusion)] <- 0
  
  grouped_data <- melt(confusion)
  setnames(grouped_data, c("genetic_assignment", "hash_assignment", "cells"))
  
  # Set assignment labels as factors so I can order them how I want
  grouped_data$genetic_assignment <- factor(grouped_data$genetic_assignment,
                                            levels = c("unassigned", "doublet",
                                                       "donor3", "donor2",
                                                       "donor1", "donor0"))
  grouped_data
}

# Load in batch 1 data
dec_data <- load_and_filter_datasets("12162021")
dec_hashing <- dec_data[[1]]
dec_genetic <- dec_data[[2]]
dec_sce <- dec_data[[3]]

# Add hash assignments to SingleCellExperiment object for UMAP plotting
# Note: for plotting we'll label multiplets as unassigned,
# but they will still be separate for the confusion matrix
colData(dec_sce) <- cbind(colData(dec_sce), "assignment_0.90" = dec_hashing$assignment_0.90)
dec_sce[, dec_sce$assignment_0.90=="Blanks" |
          dec_sce$assignment_0.90=="Multiplet"]$assignment_0.90 <- "Unassigned"
dec_hashing[dec_hashing$assignment_0.90=="Blanks", ]$assignment_0.90 <- "Unassigned"
dec_sce$assignment_0.90 <- as.factor(dec_sce$assignment_0.90)

sce_matrix <- as.data.frame(reducedDim(dec_sce, type="UMAP"))
sce_matrix$Assignment <- dec_sce$assignment_0.90
pA <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch A") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))

# Load in batch 2 data
jan_data <- load_and_filter_datasets("01132022")
jan_hashing <- jan_data[[1]]
jan_genetic <- jan_data[[2]]
jan_sce <- jan_data[[3]]

# Add hash assignments to SingleCellExperiment object for batch 2
# Note: for plotting we'll label multiplets as unassigned,
# but they will still be separate for the confusion matrix
colData(jan_sce) <- cbind(colData(jan_sce), "assignment_0.90" = jan_hashing$assignment_0.90)
jan_sce[, jan_sce$assignment_0.90=="Blanks" |
          jan_sce$assignment_0.90=="Multiplet"]$assignment_0.90 <- "Unassigned"
jan_sce$assignment_0.90 <- as.factor(jan_sce$assignment_0.90)

sce_matrix <- as.data.frame(reducedDim(jan_sce, type="UMAP"))
sce_matrix$Assignment <- jan_sce$assignment_0.90
pB <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch B") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))


# Add genetic assignments to SingleCellExperiment object for UMAP plotting
colData(dec_sce) <- cbind(colData(dec_sce), "genetic_assignment" = dec_genetic$donor_id)
dec_sce[, dec_sce$genetic_assignment=="doublet"]$genetic_assignment <- "unassigned"
dec_sce$genetic_assignment <- as.factor(dec_sce$genetic_assignment)

sce_matrix <- as.data.frame(reducedDim(dec_sce, type="UMAP"))
sce_matrix$Assignment <- dec_sce$genetic_assignment
pC <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(name = "Assignment",
                     values = c("donor0" = "#4DAC26",
                                "donor1" = "#0571B0",
                                "donor2" = "#CA0020",
                                "donor3" = "#E66101",
                                "unassigned" = "#999999"))


colData(jan_sce) <- cbind(colData(jan_sce), "genetic_assignment" = jan_genetic$donor_id)
jan_sce[, jan_sce$genetic_assignment=="doublet"]$genetic_assignment <- "unassigned"
jan_sce$genetic_assignment <- as.factor(jan_sce$genetic_assignment)

sce_matrix <- as.data.frame(reducedDim(jan_sce, type="UMAP"))
sce_matrix$Assignment <- jan_sce$genetic_assignment
pD <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(name = "Assignment",
                     values = c("donor0" = "#4DAC26",
                                "donor1" = "#0571B0",
                                "donor2" = "#CA0020",
                                "donor3" = "#E66101",
                                "unassigned" = "#999999"))


dec_confusion <- make_overlap_matrix(dec_hashing, dec_genetic, "0.90")
pE <- ggplot(dec_confusion, aes(x=hash_assignment, y=genetic_assignment, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("Hash demultiplexing assignment") + ylab("Genetic demultiplexing assignment") + 
  theme(legend.position = "None")


jan_confusion <- make_overlap_matrix(jan_hashing, jan_genetic, "0.90")
pF <- ggplot(jan_confusion, aes(x=hash_assignment, y=genetic_assignment, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("Hash demultiplexing assignment") + ylab("Genetic demultiplexing assignment") + 
  theme(legend.position = "None")

pdf("../../figures/figure2.pdf", width = 16, height = 16, family = "sans")
pA + pB + pC + pD + pE + pF +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = "A")
dev.off()



## SUPPLEMENTAL FIGURES

# Calculate assignments based on different thresholds
dec_hashing$assignment_0.85 <- set_new_threshold(dec_hashing, 0.85)
dec_hashing$assignment_0.80 <- set_new_threshold(dec_hashing, 0.80)
jan_hashing$assignment_0.85 <- set_new_threshold(jan_hashing, 0.85)
jan_hashing$assignment_0.80 <- set_new_threshold(jan_hashing, 0.80)

# Add assignment thresholds to SCE object for UMAP plotting
colData(dec_sce) <- cbind(colData(dec_sce), "assignment_0.85" = dec_hashing$assignment_0.85)
colData(dec_sce) <- cbind(colData(dec_sce), "assignment_0.80" = dec_hashing$assignment_0.80)
colData(jan_sce) <- cbind(colData(jan_sce), "assignment_0.85" = jan_hashing$assignment_0.85)
colData(jan_sce) <- cbind(colData(jan_sce), "assignment_0.80" = jan_hashing$assignment_0.80)

sce_matrix <- as.data.frame(reducedDim(dec_sce, type="UMAP"))
sce_matrix$Assignment <- dec_sce$assignment_0.85
qA <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch A, 85% Threshold") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))

sce_matrix <- as.data.frame(reducedDim(jan_sce, type="UMAP"))
sce_matrix$Assignment <- jan_sce$assignment_0.85
qB <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch B, 85% Threshold") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))

sce_matrix <- as.data.frame(reducedDim(dec_sce, type="UMAP"))
sce_matrix$Assignment <- dec_sce$assignment_0.80
qC <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch A, 80% Threshold") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))

sce_matrix <- as.data.frame(reducedDim(jan_sce, type="UMAP"))
sce_matrix$Assignment <- jan_sce$assignment_0.80
qD <- ggplot(sce_matrix, mapping=aes(x=V1, y=V2, color=Assignment)) + geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Batch B, 80% Threshold") +
  scale_color_manual(name = "Assignment",
                     values = c("Hashtag1" = "#CA0020",
                                "Hashtag2" = "#0571B0",
                                "Hashtag3" = "#4DAC26",
                                "Hashtag4" = "#E66101",
                                "Unassigned" = "#999999"))

pdf("../../figures/suppfig1.pdf", width = 15, height = 10.6, family = "sans")
qA + qB + qC + qD + plot_annotation(tag_levels = "A")
dev.off()




make_cell_type_barchart <- function(ct) {
  ct_all <- table(ct$Simplified)
  ct_all <- as.data.frame(ct_all/sum(ct_all))
  ct_all$type <- "All cells"
  
  ct_unassigned_90 <- table(ct[ct$assignment_0.90=="Unassigned",]$Simplified)
  ct_unassigned_90 <- as.data.frame(ct_unassigned_90/sum(ct_unassigned_90))
  ct_unassigned_90$type <- "Unassigned at 90% Threshold"
  
  ct_unassigned_85 <- table(ct[ct$assignment_0.85=="Unassigned",]$Simplified)
  ct_unassigned_85 <- as.data.frame(ct_unassigned_85/sum(ct_unassigned_85))
  ct_unassigned_85$type <- "Unassigned at 85% Threshold"
  
  ct_unassigned_80 <- table(ct[ct$assignment_0.80=="Unassigned",]$Simplified)
  ct_unassigned_80 <- as.data.frame(ct_unassigned_80/sum(ct_unassigned_80))
  ct_unassigned_80$type <- "Unassigned at 80% Threshold"
 
  unassigned <- rbind(ct_all, ct_unassigned_90, ct_unassigned_85, ct_unassigned_80)
  setnames(unassigned, c("cell_type", "proportion", "type"))
  
  unassigned
}

dec_ct <- load_cell_types("12162021", dec_hashing)
dec_unassigned <- make_cell_type_barchart(dec_ct)
rA <- ggplot(dec_unassigned, mapping = aes(x=type, y=proportion, fill=cell_type)) +
  geom_bar(stat = "identity") +
  ggtitle("Batch A") +
  ylab("Proportion") + 
  labs(fill = "Cell type") +
  theme(axis.title.x = element_blank())

jan_ct <- load_cell_types("01132022", jan_hashing)
jan_unassigned <- make_cell_type_barchart(jan_ct)
rB <- ggplot(jan_unassigned, mapping = aes(x=type, y=proportion, fill=cell_type)) +
  geom_bar(stat = "identity") +
  ggtitle("Batch B") +
  ylab("Proportion") + 
  labs(fill = "Cell type") +
  theme(axis.title.x = element_blank())

pdf("../../figures/suppfig2.pdf", width = 12, height = 6, family = "sans")
rA + rB + plot_annotation(tag_levels = "A")
dev.off()




get_genotype_heatmaps <- function(sample_id, type1, type2) { 
  vireo_path <- paste(data_path, "pooled_tumors", sample_id, "vireo", sep = "/")
  data1 <- fread(paste(vireo_path, type1, "donor_ids.tsv", sep = "/"))
  data2 <- fread(paste(vireo_path, type2, "donor_ids.tsv", sep = "/"))
  
  # Filter vireo to cells in sce object, all others failed miQC filtering
  sce <- readRDS(paste(local_data_path,"/sce_objects/",
                       sample_id, ".rds", sep = ""))
  data1 <- subset(data1, data1$cell %in% sce$Barcode)
  data2 <- subset(data2, data2$cell %in% sce$Barcode)
  
  setnames(data1, "donor_id", type1)
  setnames(data2, "donor_id", type2)
  
  data1 <- subset(data1, select = c("cell", type1))
  data2 <- subset(data2, select = c("cell", type2))
  
  master_list <- full_join(data1, data2)
  
  groups <- c(type1, type2)
  master_list <- subset(master_list, select = groups)
  grouped_data <- master_list %>% 
    group_by(across(all_of(groups))) %>%
    summarize(count = n()) %>% as.data.frame()
  
  # Reshape and then melt to replace NAs with 0s
  confusion <- reshape(data=grouped_data, 
                       idvar=type1,
                       v.names="count",
                       timevar=type2,
                       direction="wide")
  colnames(confusion) <- gsub("count.", "", colnames(confusion))
  confusion[is.na(confusion)] <- 0
  
  grouped_data <- melt(confusion)

  # Set assignment labels as factors so I can order them how I want
  setnames(grouped_data, c("group1", "group2", "cells"))
  grouped_data$group1 <- factor(grouped_data$group1,
                                            levels = c("donor0", "donor1",
                                                       "donor2", "donor3",
                                                       "doublet", "unassigned"))
  grouped_data$group2 <- factor(grouped_data$group2,
                                levels = c("unassigned", "doublet",
                                           "donor3", "donor2",
                                           "donor1", "donor0"))

  setnames(grouped_data, c(type1, type2, "cells"))  
  grouped_data
}

confusion1 <- get_genotype_heatmaps("12162021", "chunk_ribo", "dissociated_ribo")
sA <- ggplot(confusion1, aes(x=chunk_ribo, y=dissociated_ribo, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("chunk_ribo genotypes") + ylab("dissociated_ribo genotypes") + 
  ggtitle("Batch A") +
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5))

confusion2 <- get_genotype_heatmaps("01132022", "chunk_ribo", "dissociated_ribo")
sB <- ggplot(confusion2, aes(x=chunk_ribo, y=dissociated_ribo, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("chunk_ribo genotypes") + ylab("dissociated_ribo genotypes") + 
  ggtitle("Batch B") +
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5))

confusion3 <- get_genotype_heatmaps("12162021", "dissociated_ribo", "dissociated_polyA")
sC <- ggplot(confusion3, aes(x=dissociated_ribo, y=dissociated_polyA, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("dissociated_ribo genotypes") + ylab("dissociated_polyA genotypes") + 
  ggtitle("Batch A") +
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5))

confusion4 <- get_genotype_heatmaps("01132022", "dissociated_ribo", "dissociated_polyA")
sD <- ggplot(confusion4, aes(x=dissociated_ribo, y=dissociated_polyA, fill=cells)) +
  geom_raster() + geom_text(aes(label = cells), size = 6) +
  heatmap_scale_1d + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("dissociated_ribo genotypes") + ylab("dissociated_polyA genotypes") + 
  ggtitle("Batch B") +
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5))

pdf("../../figures/suppfig3.pdf", width = 12, height = 12, family = "sans")
sA + sB + sC + sD +
  plot_annotation(tag_levels = "A")
dev.off()
