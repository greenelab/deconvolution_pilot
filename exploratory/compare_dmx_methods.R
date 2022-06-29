# Look at results of hash and genetic demultiplexing

suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(scater)
  library(miQC)
  library(ggplot2)
  library(dplyr)
})

source("../config.R")

## Load data

# We have two sets of pooled samples, named after the date they were run (12162021 and 01132022).
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Sample ID must be provided.", call.=FALSE)
} else {
  sample_id <- args[1]
}

# Load hash and vireo (genetic) demultiplexing assignments
hashing <- fread(paste(data_path, "pooled_tumors", sample_id,
		       "Cellranger/outs/multi/multiplexing_analysis/assignment_confidence_table.csv", sep = "/"))
vireo <- fread(paste(data_path, "pooled_tumors", sample_id, "vireo/donor_ids.tsv", sep = "/"))

# Load SingleCellExperiment object for plotting
sce <- readRDS(paste("sce_objects/", sample_id, ".rds", sep = ""))

# Filter hashing and vireo to cells in sce object, all others failed miQC filtering
hashing <- subset(hashing, hashing$Barcodes %in% sce$Barcode)
vireo <- subset(vireo, vireo$cell %in% sce$Barcode)

## Explore hash demultiplexing thresholds

# Based on the hash assignment confidence table, the default cutoff for cellranger multi is to only assign a cell if
# it has >90% probability of belonging to a given hashtag. Based on a manual look at the assignment table, this seems 
# unnecessarily stringent. Here we test what number of cells we get at two different cutoff levels, 85% (moderate), 
# and 80% (liberal).

# Make new assignments based on moderate cutoffs
mod_hashtag_prob <- 0.85
other_hashtag_prob <- 0.001
moderate <- hashing
moderate[moderate$Assignment == "Unassigned" &
           moderate$`anti-human_Hashtag1` >= mod_hashtag_prob & 
           moderate$`anti-human_Hashtag2` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag3` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag1"
moderate[moderate$Assignment == "Unassigned" &
           moderate$`anti-human_Hashtag1` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag2` >= mod_hashtag_prob &
           moderate$`anti-human_Hashtag3` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag2"
moderate[moderate$Assignment == "Unassigned" &
           moderate$`anti-human_Hashtag1` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag2` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag3` >= mod_hashtag_prob &
           moderate$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag3"
moderate[moderate$Assignment == "Unassigned" &
           moderate$`anti-human_Hashtag1` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag2` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag3` <= other_hashtag_prob &
           moderate$`anti-human_Hashtag4` >= mod_hashtag_prob, ]$Assignment = "anti-human_Hashtag4"
table(moderate$Assignment)

# Make new assignments based on liberal cutoffs
lib_hashtag_prob <- 0.80
other_hashtag_prob <- 0.001
liberal <- hashing
liberal[liberal$Assignment == "Unassigned" &
          liberal$`anti-human_Hashtag1` >= lib_hashtag_prob & 
          liberal$`anti-human_Hashtag2` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag3` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag1"
liberal[liberal$Assignment == "Unassigned" &
          liberal$`anti-human_Hashtag1` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag2` >= lib_hashtag_prob &
          liberal$`anti-human_Hashtag3` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag2"
liberal[liberal$Assignment == "Unassigned" &
          liberal$`anti-human_Hashtag1` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag2` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag3` >= lib_hashtag_prob &
          liberal$`anti-human_Hashtag4` <= other_hashtag_prob, ]$Assignment = "anti-human_Hashtag3"
liberal[liberal$Assignment == "Unassigned" &
          liberal$`anti-human_Hashtag1` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag2` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag3` <= other_hashtag_prob &
          liberal$`anti-human_Hashtag4` >= lib_hashtag_prob, ]$Assignment = "anti-human_Hashtag4"
table(liberal$Assignment)


## Get overlap of demultiplexing labels

# Add assignments to sce object
assignment_cons <- hashing$Assignment
colData(sce) <- cbind(colData(sce), assignment_cons)
assignment_mod <- moderate$Assignment
colData(sce) <- cbind(colData(sce), assignment_mod)
assignment_lib <- liberal$Assignment
colData(sce) <- cbind(colData(sce), assignment_lib)
assignment_vireo <- vireo$donor_id
colData(sce) <- cbind(colData(sce), assignment_vireo)

master_list <- subset(colData(sce), select = c("Barcode", "assignment_cons", "assignment_mod",
                                               "assignment_lib", "assignment_vireo"))
master_list <- as.data.frame(master_list)

# Make pseudo-confusion matrix of cell labels
overlap_matrix <- function(hash_type) {
  groups <- c("assignment_vireo", hash_type)
  grouped_data <- master_list %>% 
    group_by(across(all_of(groups))) %>%
    summarize(count = n()) %>% as.data.frame()
  confusion <- reshape(data=grouped_data, idvar="assignment_vireo",v.names="count",timevar=hash_type,direction="wide")
  colnames(confusion) <- gsub("count.", "", colnames(confusion))
  rownames(confusion) <- confusion$assignment_vireo; confusion$assignment_vireo <- NULL
  confusion[is.na(confusion)] <- 0
  confusion <- confusion[,order(colnames(confusion))]
  
  confusion
}

overlap_matrix("assignment_cons")
overlap_matrix("assignment_mod")
overlap_matrix("assignment_lib")

## Plot assignments on UMAP

# Combine blanks, multiplets and unassigned cells into one group for clearer plotting
sce[, sce$assignment_cons=="Blanks" | sce$assignment_cons=="Multiplet"]$assignment_cons <- "Unassigned"
sce[, sce$assignment_mod=="Blanks" | sce$assignment_mod=="Multiplet"]$assignment_mod <- "Unassigned"
sce[, sce$assignment_lib=="Blanks" | sce$assignment_lib=="Multiplet"]$assignment_lib <- "Unassigned"
sce[, sce$assignment_vireo=="doublet"]$assignment_vireo <- "unassigned"

# Since one of the samples failed in one of the pools, it has no cells at more stringent thresholds and this
# messes up the indexing for color. Checking for this beforehand each time keeps unassigned cells as gray.
get_colors <- function(assignment) {
  if (length(unique(assignment)) == 5) {
    colors <- scale_color_manual(name = "assignment", values = c("#CA0020", "#0571B0", "#4DAC26", "#E66101", "#999999"))
  } else if (length(unique(assignment)) == 4) {
    colors <- scale_color_manual(name = "assignment", values = c("#0571B0", "#4DAC26", "#E66101", "#999999"))
  }
  colors
}

# Plot hash demultiplexing assignments at the different filtering thresholds
png(paste("plots/pooled/", sample_id, "_UMAP_assignment_90.png", sep=""), width = 700)
colors <- get_colors(sce$assignment_cons); plotUMAP(sce, colour_by = "assignment_cons") + colors
dev.off()

png(paste("plots/pooled/", sample_id, "_UMAP_assignment_85.png", sep=""), width = 700)
colors <- get_colors(sce$assignment_mod); plotUMAP(sce, colour_by = "assignment_mod") + colors
dev.off()

png(paste("plots/pooled/", sample_id, "_UMAP_assignment_80.png", sep=""), width = 700)
colors <- get_colors(sce$assignment_lib); plotUMAP(sce, colour_by = "assignment_lib") + colors
dev.off()

# Plot genetic demultiplexing assignments
png(paste("plots/pooled/", sample_id, "_UMAP_assignment_genetic.png", sep=""), width = 700)
colors <- get_colors(sce$assignment_vireo); plotUMAP(sce, colour_by = "assignment_vireo") + colors
dev.off()
