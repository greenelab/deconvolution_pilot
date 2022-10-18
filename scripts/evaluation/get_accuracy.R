# Not only do we have three kinds of bulk RNA-seq, we also have two kinds of
# scRNA-seq: the samples run individually and the pooled samples. Because of
# this, we can use pooled cells as our reference profiles and use the solo
# samples' cell type proportions as a kind of silver standard. We've annotated
# the cell types in the solo samples, and our assumption is that the closer the
# deconvolution results are to the single-cell proportions, the closer the
# deconvolution method gets to the true proportions.

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
cell_types <- c("NK cells", "T cells", "B cells", "Endothelial cells",
                "Macrophages", "Mast cells", "Plasma cells",
                "Fibroblasts", "Monocytes", "pDC", "DC", "ILC")

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
    x<-table(sce$cellType)/ncol(sce)
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
    if (method %in% scores){
      next
    }
    results_tmp <- fread(resultfile,header=T)
    
    results <- melt(results_tmp,id.vars = "cell_type")
    results$bulk_type <- bulk_type
    results$method <- method
    
    melted_results <- rbind(melted_results, results)
}
setnames(melted_results, "value","proportion")

# Unify cell type nomenclature across methods
melted_results[melted_results$cell_type %in% c("NK cell", "NK cells",
                                               "NK_cells", "NKcells"),]$cell_type <- "NK cells"
melted_results[melted_results$cell_type %in% c("T cell", "T cells"),]$cell_type <- "T cells"
melted_results[melted_results$cell_type %in% c("T cell CD4+", "T_cells_CD4",
                                               "CD4_Tcells"),]$cell_type <- "CD4 T cells"
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

# Filter down to cell types we have single-cell data for
melted_results <- subset(melted_results, melted_results$cell_type %in% cell_types)


# Reshape the results so each method has its own column, then add a column for single cell proportions
deconvolution <- reshape(data = melted_results,
                                 idvar = c("cell_type", "variable","bulk_type"),
                                 v.names = "proportion",
                                 timevar = "method",
                                 direction = "wide")
setnames(deconvolution, "variable", "sample")
deconvolution <- left_join(deconvolution, melted_sc)

# Subtract single cell proportions from deconvolution estimates
deconvolution$proportion.abis <- deconvolution$proportion.abis - deconvolution$proportion.sc
deconvolution$proportion.bayesprism <- deconvolution$proportion.bayesprism - deconvolution$proportion.sc
deconvolution$proportion.bisque <- deconvolution$proportion.bisque - deconvolution$proportion.sc
deconvolution$proportion.epic <- deconvolution$proportion.epic - deconvolution$proportion.sc
deconvolution$proportion.music <- deconvolution$proportion.music - deconvolution$proportion.sc
deconvolution$proportion.nnls <- deconvolution$proportion.nnls - deconvolution$proportion.sc
deconvolution$proportion.quantiseq <- deconvolution$proportion.quantiseq - deconvolution$proportion.sc
deconvolution$proportion.cibersortx <- deconvolution$proportion.cibersortx - deconvolution$proportion.sc
deconvolution$proportion.sc <- NULL

# Now that we have each proportion minus the corresponding single cell
# fraction, we melt the data again for easier plotting
remelt <- melt(deconvolution)
remelt$variable <- gsub("proportion.","",remelt$variable)

# Compare accuracy across methods, stratified by cell type and vice versa
ggplot(remelt, aes(x=cell_type, y=value, fill=variable)) + geom_boxplot() + ylab("Estimated proportion - single cell proportion") + xlab("Cell type")
ggplot(remelt, aes(x=variable, y=value, fill=cell_type)) + geom_boxplot() + ylab("Estimated proportion - single cell proportion") + xlab("Cell type")


# Make heatmap of proportion differences for each bulk type separately
chunk_ribo <- subset(remelt, remelt$bulk_type=="chunk_ribo"); chunk_ribo$bulk_type <- NULL
chunk_ribo <- as.data.frame(chunk_ribo %>% group_by(cell_type, variable) %>% summarize(value=mean(value)))
chunk_ribo <- reshape(data = chunk_ribo,
                      idvar = "cell_type",
                      timevar = "variable",
                      direction = "wide")
rownames(chunk_ribo) <- chunk_ribo$cell_type; chunk_ribo$cell_type <- NULL
setnames(chunk_ribo, c("abis","bayesprism","bisque","cibersortx","epic","music","nnls","quantiseq"))
pheatmap::pheatmap(chunk_ribo, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks=ppoints(100)*1.2-0.6)


dissociated_ribo <- subset(remelt, remelt$bulk_type=="dissociated_ribo"); dissociated_ribo$bulk_type <- NULL
dissociated_ribo <- as.data.frame(dissociated_ribo %>% group_by(cell_type, variable) %>% summarize(value=mean(value)))
dissociated_ribo <- reshape(data = dissociated_ribo,
                      idvar = "cell_type",
                      timevar = "variable",
                      direction = "wide")
rownames(dissociated_ribo) <- dissociated_ribo$cell_type; dissociated_ribo$cell_type <- NULL
setnames(dissociated_ribo, c("abis","bayesprism","bisque","cibersortx","epic","music","nnls","quantiseq"))
pheatmap::pheatmap(dissociated_ribo, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks=ppoints(100)*1.2-0.6)


dissociated_polyA <- subset(remelt, remelt$bulk_type=="dissociated_polyA"); dissociated_polyA$bulk_type <- NULL
dissociated_polyA <- as.data.frame(dissociated_polyA %>% group_by(cell_type, variable) %>% summarize(value=mean(value)))
dissociated_polyA <- reshape(data = dissociated_polyA,
                      idvar = "cell_type",
                      timevar = "variable",
                      direction = "wide")
rownames(dissociated_polyA) <- dissociated_polyA$cell_type; dissociated_polyA$cell_type <- NULL
setnames(dissociated_polyA, c("abis","bayesprism","bisque","cibersortx","epic","music","nnls","quantiseq"))
pheatmap::pheatmap(dissociated_polyA, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 18, legend = FALSE,
                   breaks=ppoints(100)*1.2-0.6)