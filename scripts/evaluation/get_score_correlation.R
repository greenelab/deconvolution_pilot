# Some deconvolution methods don't return proportions, but instead cell type
# scores that are meant to be compared across samples and don't have a real
# unit. Because of this, sometimes people doing benchmarking for these methods
# will compare the correlations between the scores and known values. For this,
# we'll look at the correlation between scores and single-cell proportions.

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
samples <- params$samples

# Which methods have values that represent true proportions (i.e. sum to 1)
proportions <- c("abis","bayesprism","bisque","cibersortx","epic","music","nnls","quantiseq")
scores <- c("consensus_tme","immucellai","mcpcounter","timer","xcell")

# Cell types of interest
cell_types <- c("NK cells", "CD4 T cells", "CD8 T cells", "Regulatory T cells", 
                "B cells", "Endothelial cells", "Macrophages", "Mast cells",
                "Plasma cells","Fibroblasts","Monocytes","pDC","DC")

# Load in labeled single cells, melted into a single data frame
melted_sc <- data.frame()
for(i in 1:length(samples)){
  sample <- samples[i]
  if(sample=="2428"){
    next
  }
  sc_file <- paste(local_data_path, "/sce_objects/", sample,
                   "_labeled.rds", sep = "")
  sce <- readRDS(sc_file)
  x<-table(sce$cellTypeGranular)/ncol(sce)
  y<-melt(x); setnames(y,"Var1","cell_type")
  y$sample <- sample
  
  melted_sc <- rbind(melted_sc, y)
}
setnames(melted_sc, "value","proportion.sc")

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
  if (method %in% proportions){
    next
  }
  results_tmp <- fread(resultfile,header=T)
  
  results <- melt(results_tmp,id.vars = "cell_type")
  results$bulk_type <- bulk_type
  results$method <- method
  
  melted_results <- rbind(melted_results, results)
}
setnames(melted_results, "value","score")
setnames(melted_results, "variable","sample")

# Unify cell type nomenclature across methods
melted_results[melted_results$cell_type %in% c("NK cell", "NK cells",
                                               "NK_cells", "NKcells", "NK"),]$cell_type <- "NK cells"
melted_results[melted_results$cell_type %in% c("T cell", "T cells"),]$cell_type <- "T cells"
melted_results[melted_results$cell_type %in% c("T cell CD4+", "T_cells_CD4",
                                               "CD4_Tcells", "CD4_T"),]$cell_type <- "CD4 T cells"
melted_results[melted_results$cell_type %in% c("CD8_Tcells", "T_cells_CD8",
                                               "CD8_T", "T cell CD8+"),]$cell_type <- "CD8 T cells"
melted_results[melted_results$cell_type %in% c("T cell regulatory (Tregs)", "T_regulatory_cells",
                                               "T regulatory cells", "nTreg"),]$cell_type <- "Regulatory T cells"
melted_results[melted_results$cell_type %in% c("B_cells", "Bcells",
                                               "B cells","B cell", "Bcell"),]$cell_type <- "B cells"
melted_results[melted_results$cell_type %in% c("Endothelial", "Endothelial cell",
                                               "Endothelial cells"),]$cell_type <- "Endothelial cells"
melted_results[melted_results$cell_type %in% c("Eosinophil", "Eosinophils"),]$cell_type <- "Eosinophils"
melted_results[melted_results$cell_type %in% c("Macrophage", "Macrophages","Macrophage/Monocyte"),]$cell_type <- "Macrophages"
melted_results[melted_results$cell_type %in% c("Mast cells", "Mast_cells",
                                               "Mast cell")]$cell_type <- "Mast cells"
melted_results[melted_results$cell_type %in% c("Plasma cells", "Plasma_cells",
                                               "B cell plasma"),]$cell_type <- "Plasma cells"
melted_results[melted_results$cell_type %in% c("Fibroblasts",
                                               "Cancer associated fibroblast"),]$cell_type <- "Fibroblasts"
melted_results[melted_results$cell_type %in% c("Monocytes", "Monocyte")]$cell_type <- "Monocytes"
melted_results[melted_results$cell_type %in% c("pDC",
                                               "Plasmacytoid dendritic cell")]$cell_type <- "pDC"
melted_results[melted_results$cell_type %in% c("Dendritic_cells","DC","Myeloid dendritic cell")]$cell_type <- "DC"

# Check how many methods report each cell type category
for(i in 1:length(table(melted_sc$cell_type))){
  celltype <- names(table(melted_sc$cell_type))[i]
  print(celltype)
  print(length(which(melted_results$cell_type==celltype))/24)
}

# Plot correlations
for(i in 1:length(scores)){
  methodname <- scores[i]
  x <- subset(melted_results, melted_results$method==methodname)
  x <- subset(x, x$sample !="2428")
  mutual_cell_types <- intersect(x$cell_type, melted_sc$cell_type)
  x <- full_join(x, melted_sc)
  x <- subset(x, x$cell_type %in% mutual_cell_types)
  
  # Mark empty categories in single cell as proportions of 0
  x[is.na(x$proportion.sc),]$proportion.sc <- 0
  
  title <- paste(methodname, ", r = ", cor(x$score, x$proportion.sc), sep = "")
  ggplot(x, mapping = aes(x=score, y=proportion.sc, color=bulk_type)) + geom_point() + ggtitle(title)
}
