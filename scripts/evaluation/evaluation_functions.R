suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(dplyr)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load in labeled single cells, melted into a single data frame
# Granular: boolean, whether to split out T cells into CD4/CD8/regulatory
load_melted_sc <- function(granular = FALSE){
    melted_sc <- data.frame()
    for(i in 1:length(samples)){
        sample <- samples[i]
        # Single cell sequencing more or less failed for this sample, with only 300
        # cells total, so we can't make reliable proportions estimates for it.
        if(sample=="2428"){
            next
        }
        sc_file <- paste(local_data_path, "/sce_objects/", sample,
                         "_labeled.rds", sep = "")
        sce <- readRDS(sc_file)
        if (granular){
            x <- table(sce$cellTypeGranular)/ncol(sce)
        } else{
            x <- table(sce$cellType)/ncol(sce)
        }
        y<-melt(x); setnames(y,"Var1","cell_type")
        y$sample <- sample
        
        melted_sc <- rbind(melted_sc, y)
    }
    setnames(melted_sc, "value","proportion.sc")
    
    melted_sc
}


# Load in pseudobulk "true fractions" (as provided by SimBu) into a dataframe
load_pseudobulk_fractions <- function(pseudobulk_types){
  melted_fractions <- data.frame()
  for(i in 1:length(pseudobulk_types)){
    pseudo_type <- pseudobulk_types[i]
    fractfile <- paste(local_data_path, "/deconvolution_input/",
                       "cell_type_fractions_", pseudo_type, ".tsv", sep = "")
    fraction <- fread(fractfile)
    setnames(fraction, "V1", "sample")
    y <- melt(fraction)
    y$bulk_type <- pseudo_type
    
    melted_fractions <- rbind(melted_fractions, y)
  }
  setnames(melted_fractions, "value", "proportion")
  setnames(melted_fractions, "variable", "cell_type")
  
  melted_fractions
}


# Load all deconvolution results into a single dataframe, with one
# row for each combination of cell type x method x bulk type x sample
load_melted_results <- function(){
    # Get the file locations of all deconvolution results
    output <- paste(local_data_path, "deconvolution_output", sep = "/")
    files <- list.files(output, full.names = T, recursive = T)
    files <- grep(".tsv", files, value = TRUE)

    melted_results <- data.frame()
    for(i in 1:length(files)){
        resultfile <- files[i]
        method <- strsplit(resultfile, split="/")
        bulk_type <- sapply(method, "[[", length(method[[1]])-1)
        method <- sapply(method, "[[", length(method[[1]]))
        method <- gsub("_results.tsv", "",  method)
        results_tmp <- fread(resultfile,header=T)
        
        results <- melt(results_tmp,id.vars = "cell_type")
        results$bulk_type <- bulk_type
        results$method <- method
        
        melted_results <- rbind(melted_results, results)
    }
    setnames(melted_results, "value","proportion")   
    
    melted_results
}


# Unify cell type nomenclature across methods
rename_cell_types <- function(results_file) {
    results_file[results_file$cell_type %in% c("NK cell", "NK cells", "NK_cells", "NKcells"),]$cell_type <- "NK cells"
    results_file[results_file$cell_type %in% c("T cell", "T cells"),]$cell_type <- "T cells"
    results_file[results_file$cell_type %in% c("T cell CD4+", "T_cells_CD4", "CD4_Tcells"),]$cell_type <- "CD4 T cells"
    results_file[results_file$cell_type %in% c("CD8_Tcells", "T_cells_CD8"),]$cell_type <- "CD8 T cells"
    results_file[results_file$cell_type %in% c("T cell regulatory (Tregs)", "T_regulatory_cells"),]$cell_type <- "T regulatory cells"
    results_file[results_file$cell_type %in% c("B_cells", "Bcells", "B cells","B cell"),]$cell_type <- "B cells"
    results_file[results_file$cell_type %in% c("Endothelial", "Endothelial cell", "Endothelial cells"),]$cell_type <- "Endothelial cells"
    results_file[results_file$cell_type %in% c("Eosinophil", "Eosinophils"),]$cell_type <- "Eosinophils"
    results_file[results_file$cell_type %in% c("Macrophage", "Macrophages"),]$cell_type <- "Macrophages"
    results_file[results_file$cell_type %in% c("Mast cells", "Mast_cells", "Mast cell")]$cell_type <- "Mast cells"
    results_file[results_file$cell_type %in% c("Plasma cells", "Plasma_cells", "B cell plasma"),]$cell_type <- "Plasma cells"
    results_file[results_file$cell_type %in% c("Fibroblasts", "Cancer associated fibroblast"),]$cell_type <- "Fibroblasts"
    results_file[results_file$cell_type %in% c("Monocytes", "Monocyte")]$cell_type <- "Monocytes"
    results_file[results_file$cell_type %in% c("pDC", "Plasmacytoid dendritic cell")]$cell_type <- "pDC"

    results_file    
}
