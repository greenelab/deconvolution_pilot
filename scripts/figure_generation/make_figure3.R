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
  library(DESeq2)
  library(rtracklayer)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load DESeq2 object
deseq_path <- paste(local_data_path, "deseq2_output", sep = "/")
dds <- readRDS(paste(deseq_path, "chunk_vs_dissociated_data.rds", sep = "/"))

# Rlog transformation to adjust for heteroskedasticity
rld <- rlog(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(rld, intgroup = c("condition", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pA <- ggplot(pcaData, aes(PC1, PC2, color = sample, shape = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()

# Get DESeq2 results
res <- results(dds)

# Filter down to protein coding genes
genefile <- paste(data_path,"index/refdata-gex-GRCh38-2020-A/genes/genes.gtf",sep = "/")
gff <- readGFF(genefile)
protein_coding <- subset(gff, gff$gene_type=="protein_coding")
res <- subset(res, rownames(res) %in% protein_coding$gene_name)

# Set genes with p value of 0 to 300
res_df <- as.data.frame(res)
res_df[is.na(res_df$padj),]$padj <- 1
res_df$group <- "Other"

# Get erythrocyte genes from Xie et al 2020 figure S6
rbc_genes <- fread(paste(local_data_path, "miscellaneous",
                         "erythrocyte_genes.tsv", sep = "/"),
                   header = FALSE)
rbc_genes <- rbc_genes$V1
res_df[rownames(res_df) %in% rbc_genes,]$group <- "RBCs"

# Get marker gene sets from Emont et al 2022 table S1
emont_genes <- fread(paste(local_data_path, "miscellaneous", "emont_genes.tsv", sep = "/"))

# Get adipocyte genes
adipo_genes <- subset(emont_genes, emont_genes$cluster=="adipocyte" & emont_genes$avg_log2FC>1)$gene
res_df[rownames(res_df) %in% adipo_genes,]$group <- "Adipocytes"

# Get endo genes
endo_genes <- subset(emont_genes, emont_genes$cluster=="endothelial" & emont_genes$avg_log2FC>1)$gene
res_df[rownames(res_df) %in% endo_genes,]$group <- "Endothelial cells"

# Get macrophage genes
#macro_genes <- subset(emont_genes, emont_genes$cluster=="macrophage" & emont_genes$avg_log2FC>1)$gene
#res_df[rownames(res_df) %in% macro_genes,]$group <- "Macrophages"


# Rearrange so interesting genes are plotted on top
res_df$group <- factor(res_df$group, levels = c("Other", "Endothelial cells", "Adipocytes", "RBCs"))
res_df <- res_df[order(res_df$group), ]

# Make volcano plot
pB <- ggplot(res_df, mapping = aes(x = log2FoldChange,
                                   y = -log10(padj),
                                   color = group)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  scale_color_manual(name = "Gene Type",
                     values = c("Adipocytes" = "#7CAE00",
                                "RBCs" = "#00BFC4",
                                "Endothelial cells" = "#C77CFF",
#                                "Macrophages" = "red",
                                "Other" = "#999999"))


pdf("../../figures/figure3.pdf", width = 12, height = 8, family = "sans")
pA + pB + pC + pD +
  plot_annotation(tag_levels = "A")
dev.off()
