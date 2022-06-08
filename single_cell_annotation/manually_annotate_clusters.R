# This is the code I

suppressPackageStartupMessages({
  library(data.table)
  library(scater)
  library(scran)
  library(batchelor)
  library(igraph)
})

pool1 <- readRDS("sce_objects/12162021.rds")
pool1$Pool <- "12162021"
pool2 <- readRDS("sce_objects/01132022.rds")
pool2$Pool <- "01132022"

sce <- cbind(pool1, pool2)

set.seed(531)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)
plotUMAP(sce, colour_by = "Pool")

set.seed(531)
mnn <- fastMNN(sce,  batch = sce$Pool,
               BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam())

reducedDim(sce,  "MNN") <- reducedDim(mnn,  "corrected")
rm(mnn)

set.seed(531)
sce <- runUMAP(sce,
               dimred = "MNN",
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)
plotUMAP(sce, colour_by = "Pool")

plotUMAP(sce, colour_by = "PTPRC")
plotUMAP(sce, colour_by = "EPCAM")
plotUMAP(sce, colour_by = "PAX8")
plotUMAP(sce, colour_by = "ACTA2")

g <- buildSNNGraph(sce,
                   k = 50,   # higher = bigger clusters
                   BNPARAM = BiocNeighbors::AnnoyParam(),
                   BPPARAM = BiocParallel::MulticoreParam())
clusters <- as.factor(igraph::cluster_louvain(g)$membership)
sce$clusters <- clusters
saveRDS(sce,  file = "sce_objects/pooled_clustered_50.rds")

plotUMAP(sce, colour_by = "clusters")


getOnlySignificantGenes <- function(mylist) {
  for (i in 1:length(mylist)) {
    x <- mylist[[i]]
    x <- subset(x, select = c("p.value", "FDR"))
    y <- subset(x, x$FDR <= 0.05)
    if (nrow(y) > 0) {
      y <- y[order(y$FDR), ]
    }
    mylist[[i]] <- y
  }
  as.list(mylist)
}

# Let"s start by seeing if we can pin down any particular clusters with global

markers <- findMarkers(sce,  groups = sce$clusters,
                       block = sce$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.5,
                       pval.type = "all",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

markers <- findMarkers(sce,  groups = sce$clusters,
                       block = sce$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.3,
                       pval.type = "all",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

reference_genes <- rownames(sce)

markers <- findMarkers(sce,  groups = sce$clusters,
                       block = sce$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.5,
                       pval.type = "some",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

# Try webGestaltR
for (i in 1:length(unique(sce$clusters))) {
  cluster_genes <- sig_genes[[i]]
  cluster_pathways <- WebGestaltR(enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                  interestGene = rownames(cluster_genes),
                                  interestGeneType = "genesymbol",
                                  referenceGene = reference_genes,
                                  referenceGeneType = "genesymbol",
                                  isOutput = FALSE)
  cluster_pathways$userId <- NULL; cluster_pathways$overlapId <- NULL
}

# After some grueling googling, I have a tentative identification for most
# markers. I"m going to try breaking them up by these definitions and seeing
# what happens

fibroblasts <- sce[, sce$clusters %in% c(2, 6, 10)]
g <- buildSNNGraph(fibroblasts,
                   k = 50,   # higher = bigger clusters
                   BNPARAM = BiocNeighbors::AnnoyParam(),
                   BPPARAM = BiocParallel::MulticoreParam())
clusters <- as.factor(igraph::cluster_louvain(g)$membership)
fibroblasts$clusters <- clusters

plotUMAP(fibroblasts, colour_by = "clusters")

markers <- findMarkers(fibroblasts,  groups = fibroblasts$clusters,
                       block = fibroblasts$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.5,
                       pval.type = "some",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

# Looking at the somewhat enigmatic immune clusters (5 was clearly macrophages)
immune <- sce[, sce$clusters %in% c(1, 3, 8, 11, 12)]
g <- buildSNNGraph(immune,
                   k = 50,   # higher = bigger clusters
                   BNPARAM = BiocNeighbors::AnnoyParam(),
                   BPPARAM = BiocParallel::MulticoreParam())
clusters <- as.factor(igraph::cluster_louvain(g)$membership)
immune$clusters <- clusters

plotUMAP(immune, colour_by = "clusters")

markers <- findMarkers(immune,  groups = immune$clusters,
                       block = immune$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.3,
                       pval.type = "all",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

# We have to go deeper, looking at the two clusters that were cluster 8 in the main
apcs <- immune[, immune$clusters %in% c(5, 8)]
markers <- findMarkers(apcs,  groups = apcs$clusters,
                       block = apcs$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.5,
                       pval.type = "all",
                       BPPARAM = BiocParallel::MulticoreParam())

sig_genes <- getOnlySignificantGenes(markers)
sig_genes

# Tcells maybe?
tcells <- sce[, sce$clusters %in% c(1, 3, 12)]
markers <- findMarkers(tcells,  groups = tcells$clusters,
                       block = tcells$Pool,  # use to get within-donor DE
                       direction = "up",  lfc = 1.5,
                       pval.type = "any",
                       BPPARAM = BiocParallel::MulticoreParam())
sig_genes <- getOnlySignificantGenes(markers)
sig_genes
