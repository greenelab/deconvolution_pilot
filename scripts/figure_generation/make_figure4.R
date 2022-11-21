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
figure_path <- params$figure_path
samples <- params$samples

source("figure_utils.R")

# Load DESeq2 object
deseq_path <- paste(local_data_path, "deseq2_output", sep = "/")
dds <- readRDS(paste(deseq_path, "ribo_vs_polyA_data.rds", sep = "/"))

# Rlog transformation to adjust for heteroskedasticity
rld <- rlog(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(rld, intgroup = c("condition", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pA <- ggplot(pcaData, aes(PC1, PC2, color = sample, shape = condition)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  labs(color = "Sample", shape = "Library prep") + 
  coord_fixed()  +
  scale_color_manual(values = colors_samples)

# Get DESeq2 results
res <- results(dds)

# Filter down to protein coding genes
genefile <- paste(data_path,"index/refdata-gex-GRCh38-2020-A/genes/genes.gtf",sep = "/")
gff <- readGFF(genefile)
protein_coding <- subset(gff, gff$gene_type=="protein_coding")
res <- subset(res, rownames(res) %in% protein_coding$gene_name)

# Set genes with p value of 0 to 300
res_df <- as.data.frame(res)
res_df[res_df$padj==0,]$padj <- 1e-300
res_df$group <- "Other"

# Get poly-A negative genes from Yang et al 2011
polyA_neg <- fread(paste(local_data_path, "miscellaneous",
                         "polyA_negative_genes.tsv", sep = "/"),
                   header = FALSE)
polyA_neg <- polyA_neg$V1
res_df[rownames(res_df) %in% polyA_neg,]$group <- "Other polyA(-)"

# Split out histone genes
hist_genes <- grep("HIST", rownames(dds), value = T)
res_df[rownames(res_df) %in% hist_genes,]$group <- "Histones"

# Get mitochondrial genes
mt_genes <- grep("MT-", rownames(dds), value = T)
res_df[rownames(res_df) %in% mt_genes,]$group <- "MT Genes"

# Rearrange so interesting genes are plotted on top
res_df$group <- factor(res_df$group, levels = c("Other","Histones", "Other polyA(-)", "MT Genes"))
res_df <- res_df[order(res_df$group), ]

# Make volcano plot
pB <- ggplot(res_df, mapping = aes(x = log2FoldChange,
                                   y = -log10(padj),
                                   color = group)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  scale_color_manual(name = "Gene set", values = colors_genesets,
                     limits = c("Histones", "Other polyA(-)", "MT Genes", "Other")) +
  xlab("log2 fold change") + ylab("-log10 adjusted p-value")

# Get sums of all histone genes
hist_expr <- as.data.frame(colSums(assay(dds[hist_genes,])))
colnames(hist_expr) <- "counts"
hist_expr$id <- colData(dds)$id
hist_expr$sample <- gsub("_.*","", hist_expr$id)
hist_expr$condition <- gsub(".*_", "", hist_expr$id)

pC <- ggplot(hist_expr, aes(x=condition, y=counts, group=sample, color=sample)) +
  geom_point() +
  scale_y_log10() +
  geom_line() +
  theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Status", y = "Read counts", color = "Sample") +
  ggtitle("Histone genes") +
  scale_color_manual(values = colors_samples) +
  annotation_logticks(sides = "l")


# Get sums of all mitochondrial genes
mito_expr <- as.data.frame(colSums(assay(dds[mt_genes,])))
colnames(mito_expr) <- "counts"
mito_expr$id <- colData(dds)$id
mito_expr$sample <- gsub("_.*","", mito_expr$id)
mito_expr$condition <- gsub(".*_", "", mito_expr$id)

pD <- ggplot(mito_expr, aes(x=condition, y=counts, group=sample, color=sample)) +
  geom_point() + 
  scale_y_log10() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Status", y = "Read counts", color = "Sample") +
  ggtitle("Mitochondrial genes") +
  scale_color_manual(values = colors_samples) +
  annotation_logticks(sides = "l")


pdf(paste(figure_path, "figure4.pdf", sep = "/"), width = 16, height = 10.67, family = "sans")
pA + pB + pC + pD +
  plot_layout(nrow = 2, heights = c(1,1)) +
  plot_annotation(tag_levels = "A")
dev.off()
