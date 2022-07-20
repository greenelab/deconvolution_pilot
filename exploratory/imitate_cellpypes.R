# After showing the less-than-exciting initial deconvolution results, Stephanie
# suggested the cell type signal might be getting muddled, so that if we
# filtered down to the cells we're most certain are representative of that cell
# type, we'd get better deconvolution results. Here I'm going to try to do that.
# This is loosely based on cellpypes https://github.com/FelixTheStudent/cellpypes,
# which is structured around Seurat objects and was giving me more trouble than
# it was worth, so I have reconstructed the basic logic manually.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(SingleCellExperiment)
  library(scater)
})

# Load sce data
sce <- readRDS("../sce_objects/pooled_clustered_50.rds")

# Rerun UMAP without the custom graphical parameters
# to have more separation between clusters
set.seed(531)
sce <- runUMAP(sce,
               dimred = "MNN",
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam())

# Split out markers of cell type. This was iterated over countless times
# and is not perfect, but here's the gating used for the results I first
# shared with Casey and Stephanie on 6/27/2022.
sce$epi_markers <- ifelse(assay(sce, "counts")["SLPI", ] > 2 &
                            assay(sce, "counts")["WFDC2", ] > 4 &
                            assay(sce, "counts")["ACTA2", ] < 1 &
                            assay(sce, "counts")["VWF", ] < 1 &
                            assay(sce, "counts")["PTPRC", ] < 1 &
                            reducedDim(sce, type = "UMAP")[, 2] < -5, 1, 0)
epi_selected <- sce[, sce$epi_markers == 1]

sce$plasma_markers <- ifelse(assay(sce, "counts")["JCHAIN", ] > 1 &
                               assay(sce, "counts")["MZB1", ] > 1 &
                               assay(sce, "counts")["PTPRC", ] < 1 &
                               assay(sce, "counts")["VWF", ] < 1 &
                               reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
plasma_selected <- sce[, sce$plasma_markers == 1]

sce$endo_markers <- ifelse(assay(sce, "counts")["VWF", ] > 1 &
                             assay(sce, "counts")["GNG11", ] > 0 &
                             assay(sce, "counts")["ACTA2", ] < 1 &
                             assay(sce, "counts")["PTPRC", ] < 1 &
                             reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
endo_selected <- sce[, sce$endo_markers == 1]

sce$bcell_markers <- ifelse(assay(sce, "counts")["MS4A1", ] > 1 &
                              assay(sce, "counts")["CD3E", ] < 1 &
                              assay(sce, "counts")["CD3D", ] < 1 &
                              reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
bcell_selected <- sce[, sce$bcell_markers == 1]

sce$pdc_markers <- ifelse(assay(sce, "counts")["PLD4", ] > 2 &
                            assay(sce, "counts")["CD3E", ] < 1 &
                            assay(sce, "counts")["MS4A1", ] < 1 &
                            assay(sce, "counts")["LYZ", ] < 1 &
                            assay(sce, "counts")["AIF1", ] < 1 &
                            assay(sce, "counts")["IL1B", ] < 1 &
                            reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
pdc_selected <- sce[, sce$pdc_markers == 1]

sce$macro_markers <- ifelse(assay(sce, "counts")["LYZ", ] > 1 &
                              assay(sce, "counts")["AIF1", ] > 1 &
                              assay(sce, "counts")["CD3E", ] < 1 &
                              assay(sce, "counts")["CD3D", ] < 1 &
                              assay(sce, "counts")["EPCAM", ] < 1 &
                              assay(sce, "counts")["VWF", ] < 2 &
                              assay(sce, "counts")["JCHAIN", ] < 2 &
                              reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
macro_selected <- sce[, sce$macro_markers == 1]

sce$tcell_markers <- ifelse(assay(sce, "counts")["CD3E", ] > 0 &
                              assay(sce, "counts")["CD3D", ] > 0 &
                              assay(sce, "counts")["CD2", ] > 0 &
                              assay(sce, "counts")["EPCAM", ] < 1 &
                              assay(sce, "counts")["CD68", ] < 1 &
                              assay(sce, "counts")["MS4A1", ] < 1 &
                              reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
tcell_selected <- sce[, sce$tcell_markers == 1]

sce$nk_markers <- ifelse(assay(sce, "counts")["NKG7", ] > 1 &
                           assay(sce, "counts")["GNLY", ] > 1 &
                           assay(sce, "counts")["CD3E", ] < 1 &
                           assay(sce, "counts")["CD3D", ] < 1 &
                           assay(sce, "counts")["CD4", ] < 1 &
                           assay(sce, "counts")["CD8A", ] < 1 &
                           reducedDim(sce, type = "UMAP")[, 1] < 6, 1, 0)
nk_selected <- sce[, sce$nk_markers == 1]

sce$fibro_markers <- ifelse(assay(sce, "counts")["DCN", ] > 1 &
                              assay(sce, "counts")["PTPRC", ] < 1 &
                              assay(sce, "counts")["CD3E", ] < 1 &
                              assay(sce, "counts")["CD3D", ] < 1 &
                              assay(sce, "counts")["CD4", ] < 1 &
                              assay(sce, "counts")["CD8A", ] < 1 &
                              assay(sce, "counts")["WFDC2", ] < 3 &
                              assay(sce, "counts")["VWF", ] < 1 &
                              assay(sce, "counts")["LYZ", ] < 2 &
                              reducedDim(sce, type = "UMAP")[, 1] > 6, 1, 0)
fibro_selected <- sce[, sce$fibro_markers == 1]

# Check if any cells are assigned to >1 cell type based on this gating
sce$count <- sce$endo_markers +
  sce$epi_markers +
  sce$plasma_markers +
  sce$tcell_markers +
  sce$nk_markers +
  sce$pdc_markers +
  sce$bcell_markers +
  sce$macro_markers +
  sce$fibro_markers
table(sce$count)

# Assign cell type labels
sce$cell_type <- "Unassigned"
sce[, sce$endo_markers == 1 & sce$count == 1]$cell_type <- "Endothelial"
sce[, sce$epi_markers == 1 & sce$count == 1]$cell_type <- "Epithelial"
sce[, sce$plasma_markers == 1 & sce$count == 1]$cell_type <- "Plasma cells"
sce[, sce$tcell_markers == 1 & sce$count == 1]$cell_type <- "T cells"
sce[, sce$nk_markers == 1 & sce$count == 1]$cell_type <- "NK cells"
sce[, sce$pdc_markers == 1 & sce$count == 1]$cell_type <- "pDC"
sce[, sce$bcell_markers == 1 & sce$count == 1]$cell_type <- "B cells"
sce[, sce$macro_markers == 1 & sce$count == 1]$cell_type <- "Macrophages"
sce[, sce$fibro_markers == 1 & sce$count == 1]$cell_type <- "Fibroblasts"

# Generate heatmap using markers that we used and expect to be
# indicative of only one cell type (e.g. not CD45/PTPRC).
markers_used <- c("DCN", "CD3E", "CD3D", "WFDC2", "VWF", "LYZ",
                  "GNLY", "NKG7", "CD2", "CD68", "MS4A1",
                  "JCHAIN", "AIF1", "IL1B", "ACTA2", "GNG11",
                  "SLPI", "MZB1")
plotHeatmap(assigned, features = markers_used,
            color_columns_by = "cell_type", zlim = c(0, 6))

# Save object with only assigned cells
assigned <- sce[, sce$cell_type != "Unassigned"]
saveRDS(assigned, "../sce_objects/manually_assigned.rds")
