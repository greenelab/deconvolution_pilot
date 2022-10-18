# BayesPrism (10.1038/s43018-022-00356-3)
# https://github.com/Danko-Lab/BayesPrism

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(BayesPrism)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Get bulk counts matrix
# Note: data_path is loaded from config.R
bulk_matrix <- matrix()
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id, bulk_type,
                         "STAR/ReadsPerGene.out.tab", sep = "/")
  bulk_tmp <- fread(bulk_tmp_file)
  bulk_tmp <- bulk_tmp[grep("ENSG", bulk_tmp$V1), ]
  setnames(bulk_tmp, "V2", sample_id)
  if (i == 1) {
    bulk_matrix <- as.matrix(bulk_tmp[, 2])
    rownames(bulk_matrix) <- bulk_tmp$V1
  } else {
    bulk_matrix <- cbind(bulk_matrix, as.matrix(bulk_tmp[, 2]))
  }
}
bulk_matrix <- t(bulk_matrix)
rm(i, sample_id, bulk_tmp_file, bulk_tmp); gc()

# Get single cell counts matrix
# Note: local_data_path is loaded from config.R
sce <- readRDS(paste(local_data_path, "deconvolution_input", "labeled_cell_state_profile.rds", sep = "/"))
rownames(sce) <- rowData(sce)$ID
colnames(sce) <- sce$unique_barcode
single_cell_matrix <- t(as.matrix(assay(sce)))

# Get cell type and cell state labels
cell_types <- sce$cellType
cell_states <- sce$cellState

# Cleanup genes
sc.dat.filtered <- cleanup.genes(input = single_cell_matrix, input.type = "count.matrix",
                                 species="hs", gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1"), exp.cells = 5)

# Filter to protein coding genes for faster computation
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")

# Make prism object
myPrism <- new.prism(reference = sc.dat.filtered.pc,
                     mixture = bulk_matrix,
                     input.type = "count.matrix",
                     cell.type.labels = cell_types,
                     cell.state.labels = cell_states,
                     key="Epithelial cells")

# Run BayesPrism
bp.res <- run.prism(prism = myPrism, n.cores=6)
theta <- get.fraction(bp = bp.res,
                      which.theta = "final",
                      state.or.type = "type")

# Save BayesPrism object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "bayesprism_results_full.rds", sep = "/")
saveRDS(bp.res, file = object_file)

# Format text version of proportion estimates
theta <- as.data.frame(t(theta))
theta <- cbind(rownames(theta), theta)
colnames(theta) <- c("cell_type", samples)

# Save proportion estimates
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "bayesprism_results.tsv", sep = "/")
write.table(theta, file = text_file, sep = "\t", quote = F, row.names = F)
