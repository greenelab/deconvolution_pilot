# Most deconvolution methods that employ single-cell data require a reference
# profile with each cell labeled by cell type. This script generates that,
# using the intersection of cellTypist and unsupervised clustering results
# to determine cell labels. Note: this is using data from the pooled and
# demultiplexed runs, meaing that we can use the individually sequenced single
# cell runs to generate pseudobulk files to evaluate without testing on our
# training data.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load single cell data
sce <- readRDS(paste(local_data_path, 
                     "sce_objects/pooled_clustered_50.rds", sep = "/"))
sce$unique_barcode <- paste(sce$Pool, sce$Barcode, sep = "-")

# Load celltypist results
ct <- fread(paste(local_data_path,
                  "celltypist_output/pooled_predicted_labels.csv", sep = "/"))


# Consolidate groups for easier deconvolution
label_table <- fread(paste(local_data_path,
                     "celltypist_output/simplified_labels.tsv", sep = "/"))
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)

# Get overlap of cellTypist and manual annotation of unsupervised clustering
ct$keep <- FALSE
ct$clusters <- sce$clusters
ct[ct$clusters %in% c(1, 3, 12) & ct$Simplified %in% c("T cells", "NK cells", "ILC", "Mast cells"), ]$keep <- TRUE
ct[ct$clusters %in% c(2, 6, 10) & ct$Simplified == "Fibroblasts", ]$keep <- TRUE
ct[ct$clusters == 4 & ct$Simplified == "Endothelial cells", ]$keep <- TRUE
ct[ct$clusters == 5 & ct$Simplified %in% c("Monocytes", "Macrophages", "DC"), ]$keep <- TRUE
ct[ct$clusters %in% c(7, 9) & ct$Simplified == "Epithelial cells", ]$keep <- TRUE
ct[ct$clusters == 8 & ct$Simplified %in% c("pDC", "B cells"), ]$keep <- TRUE
ct[ct$clusters == 11 & ct$Simplified == "Plasma cells", ]$keep <- TRUE

# Keep cells within overlap
sce$cellType <- ct$Simplified
sce <- sce[, ct$keep]

## Genetic demultiplexing

# Load in cell labels from genetic demultiplexing
labels_dec <- fread(paste(data_path, "pooled_tumors/12162021/vireo/donor_ids.tsv", sep = "/"))
labels_dec$unique_barcode <- paste("12162021", labels_dec$cell, sep = "-")
labels_jan <- fread(paste(data_path, "pooled_tumors/01132022/vireo/donor_ids.tsv", sep = "/"))
labels_jan$unique_barcode <- paste("01132022", labels_jan$cell, sep = "-")
labels <- rbind(labels_dec, labels_jan)

# Remove doublets and unassigned cells
sce_genetic <- sce
labels <- labels[labels$unique_barcode %in% sce$unique_barcode, ]
sce_genetic$donor_id <- labels$donor_id
sce_genetic <- sce_genetic[, sce_genetic$donor_id %in% c("donor0", "donor1", "donor2", "donor3")]

# Save object
outfile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
saveRDS(sce_genetic, file = outfile)

## Hash demultiplexing

# Load hashing demultiplexing assignment info to get barcodes to keep
dec_hashing <- fread(paste(data_path, "pooled_tumors", "12162021",
                       "Cellranger/outs/multi/multiplexing_analysis",
                       "assignment_confidence_table.csv", sep = "/"))
dec_hashing$unique_barcode <- paste("12162021", dec_hashing$Barcodes, sep = "-")
dec_assigned <- subset(dec_hashing, dec_hashing$Assignment == "anti-human_Hashtag1" |
                         dec_hashing$Assignment == "anti-human_Hashtag2" |
                         dec_hashing$Assignment == "anti-human_Hashtag3" |
                         dec_hashing$Assignment == "anti-human_Hashtag4")$unique_barcode

jan_hashing <- fread(paste(data_path, "pooled_tumors", "01132022",
                           "Cellranger/outs/multi/multiplexing_analysis",
                           "assignment_confidence_table.csv", sep = "/"))
jan_hashing$unique_barcode <- paste("01132022", jan_hashing$Barcodes, sep = "-")
jan_assigned <- subset(jan_hashing, jan_hashing$Assignment == "anti-human_Hashtag1" |
                         jan_hashing$Assignment == "anti-human_Hashtag2" |
                         jan_hashing$Assignment == "anti-human_Hashtag3" |
                         jan_hashing$Assignment == "anti-human_Hashtag4")$unique_barcode

assigned <- c(dec_assigned, jan_assigned)

# Remove doublets and unassigned cells
sce_hash <- sce
sce_hash <- sce_hash[, sce$unique_barcode %in% assigned]

# Save object
outfile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile_hashing.rds", sep = "/")
saveRDS(sce_hash, file = outfile)


## Stress test of extreme downsampling of cells for reference profile

# Calculate proportion
table(sce$cellType)

proportions <- table(sce$cellType)/length(sce$cellType)
print(proportions)

cell_counts <- c(2000, 1000, 500, 200)

for (i in cell_counts) {
  proportions_sim <- round(proportions * i)
	proportions_sim[proportions_sim < 1] <- 1

	print(proportions_sim)

	# loop over all wanted cell-types and sample to have the final amount
	sampled_cells <- lapply(seq_along(proportions_sim), function(x) { 
	  # get all cells with the current type
		cells_of_type_x <- data.frame(colData(sce)[colData(sce)[["cellType"]] == names(proportions_sim[x]), ])
		if (proportions_sim[x] == 0) {
			cells <- c()
		} else {
			# how many cells of this type do we need?
			cells <- dplyr::slice_sample(cells_of_type_x, n = proportions_sim[x], replace = FALSE)
			cells <- cells[["Barcode"]]
		}
		return(cells)
	})
	sce_sampled <- sce[, sce$Barcode %in% unlist(sampled_cells)]

	outfile <- paste(local_data_path, "/deconvolution_input/labeled_single_cell_profile_sim", i, ".rds", sep = "")
	saveRDS(sce_sampled, file = outfile)
}


# Save text version of object for scanpy (leaving in for if I add more methods)
#rownames(sce) <- rowData(sce)$ID
#counts <- as.matrix(assay(sce))
#colnames(counts) <- sce$unique_barcode
#textfile <- paste(local_data_path, "deconvolution_input", "single_cell_matrix.tsv", sep = "/")
#write.table(counts, textfile, sep = "\t", quote = F)

#labelfile <- paste(local_data_path, "deconvolution_input", "cell_labels.tsv", sep = "/")
#write(sce$cellType, labelfile, sep = "\n")
