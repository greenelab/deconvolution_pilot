# Since we're running three kinds of bulk RNA-seq (one with tumor chunks and
# ribo depletion, one with dissociated cells and ribo depletion, one with
# dissociated cells and poly-A capture), we have three sets of deconvolution
# results for each tumor. We want to see how much proportion estimates vary
# based on the type of bulk data. Our assumption is that if a method gives
# results that are more similar/have a lower variance, they will be more
# robust to technical noise.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("abis","bayesprism","bisque","cibersortx","epic","music","nnls","quantiseq")
scores <- c("consensus_tme","immucellai","mcpcounter","timer","xcell")

# Cell types of interest
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells",
                "Macrophages", "Mast cells", "Plasma cells",
                "Fibroblasts", "Monocytes", "pDC", "DC")


# Get the file locations of all deconvolution results
output <- paste(local_data_path, "deconvolution_output", sep="/")
files <- list.files(output, full.names = T, recursive = T)
files <- grep(".tsv", files, value=TRUE)

# Load all deconvolution results into a single dataframe, with one
# row for each combination of cell type x method x bulk type x sample
melted_results <- data.frame()
for(i in 1:length(files)){
  resultfile <- files[i]
  method <- strsplit(resultfile, split="/")
  bulk_type <- sapply(method, "[[", length(method[[1]])-1)
  method <- sapply(method, "[[", length(method[[1]]))
  method <- gsub("_results.tsv", "",  method)
  # If the method's results aren't constrained to (0, 1) and/or the total
  # summing to 1, their variance won't be comparable. We'll exclude them.
  if (method %in% scores){
    next
  }
  results_tmp <- fread(resultfile, header=T)
  
  results <- melt(results_tmp,id.vars = "cell_type")
  results$bulk_type <- bulk_type
  results$method <- method
  
  melted_results <- rbind(melted_results, results)
}

# Unify cell type nomenclature across methods
melted_results[melted_results$cell_type %in% c("NK cell", "NK cells",
                                               "NK_cells", "NKcells"), ]$cell_type <- "NK cells"
melted_results[melted_results$cell_type %in% c("T cell", "T cells"),]$cell_type <- "T cells"
melted_results[melted_results$cell_type %in% c("T cell CD4+", "T_cells_CD4",
                                               "CD4_Tcells"), ]$cell_type <- "CD4 T cells"
melted_results[melted_results$cell_type %in% c("CD8_Tcells", "T_cells_CD8"),]$cell_type <- "CD8 T cells"
melted_results[melted_results$cell_type %in% c("T cell regulatory (Tregs)",
                                               "T_regulatory_cells"),]$cell_type <- "T regulatory cells"
melted_results[melted_results$cell_type %in% c("B_cells", "Bcells",
                                               "B cells","B cell"),]$cell_type <- "B cells"
melted_results[melted_results$cell_type %in% c("Endothelial", "Endothelial cell",
                                               "Endothelial cells"),]$cell_type <- "Endothelial cells"
melted_results[melted_results$cell_type %in% c("Eosinophil", "Eosinophils"),]$cell_type <- "Eosinophils"
melted_results[melted_results$cell_type %in% c("Macrophage", "Macrophages"),]$cell_type <- "Macrophages"
melted_results[melted_results$cell_type %in% c("Mast cells", "Mast_cells",
                                               "Mast cell")]$cell_type <- "Mast cells"
melted_results[melted_results$cell_type %in% c("Plasma cells", "Plasma_cells",
                                               "B cell plasma"),]$cell_type <- "Plasma cells"
melted_results[melted_results$cell_type %in% c("Fibroblasts",
                                               "Cancer associated fibroblast"),]$cell_type <- "Fibroblasts"
melted_results[melted_results$cell_type %in% c("Monocytes", "Monocyte")]$cell_type <- "Monocytes"
melted_results[melted_results$cell_type %in% c("pDC",
                                               "Plasmacytoid dendritic cell")]$cell_type <- "pDC"

# For each method x cell type x sample combination, calculate variance across bulk types
results <- melted_results %>% group_by(method,cell_type,variable) %>% summarize(variance=var(value))
results <- subset(results, results$cell_type %in% cell_types)

# Compare variance across methods overall
ggplot(results, mapping = aes(x=method, y=variance)) + geom_boxplot()

# Compare variance across methods, stratified by cell type and vice versa
ggplot(results, mapping = aes(x=method, y=log(variance), fill=cell_type)) + geom_boxplot()
ggplot(results, mapping = aes(x=cell_type, y=log(variance), fill=method)) + geom_boxplot()