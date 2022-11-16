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
  library(WebGestaltR)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

source("figure_utils.R")

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
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  labs(color = "Sample", shape = "Status") + 
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
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  scale_color_manual(name = "Gene set", values = colors_genesets, 
                     limits = c("Adipocytes", "RBCs", "Endothelial cells", "Other")) +
  xlab("log2 fold change") + ylab("-log10 adjusted p-value")

# Plot hemoglobin genes
hemo_genes <- c("HBA1","HBA2","HBB")
d1 <- plotCounts(dds, gene=hemo_genes[1], intgroup=c("condition","sample"), returnData=TRUE)
d1$Gene <- hemo_genes[1]
d2 <- plotCounts(dds, gene=hemo_genes[2], intgroup=c("condition","sample"), returnData=TRUE)
d2$Gene <- hemo_genes[2]
d3 <- plotCounts(dds, gene=hemo_genes[3], intgroup=c("condition","sample"), returnData=TRUE)
d3$Gene <- hemo_genes[3]
d <- rbind(d1, d2, d3)

pC <- ggplot(d, aes(x=condition, y=count, group=sample, color=sample)) +
    geom_point() + scale_y_log10() + geom_line() +
    facet_wrap(~Gene) +
    theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 0.5)) +
    labs(x = "Status", y = "Read counts", color = "Sample") +
    scale_color_manual(values = colors_samples)

  
# Plot adipocyte genes
selected_adipo_genes <- c("ADIPOQ", "CIDEC", "PLIN1")
d1 <- plotCounts(dds, gene=selected_adipo_genes[1], intgroup=c("condition","sample"), returnData=TRUE)
d1$Gene <- selected_adipo_genes[1]
d2 <- plotCounts(dds, gene=selected_adipo_genes[2], intgroup=c("condition","sample"), returnData=TRUE)
d2$Gene <- selected_adipo_genes[2]
d3 <- plotCounts(dds, gene=selected_adipo_genes[3], intgroup=c("condition","sample"), returnData=TRUE)
d3$Gene <- selected_adipo_genes[3]
d <- rbind(d1, d2, d3)

pD <- ggplot(d, aes(x=condition, y=count, group=sample, color=sample)) +
    geom_point() + scale_y_log10() + geom_line() +
    facet_wrap(~Gene) +
    theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Status", y = "Read counts", color = "Sample") +
    scale_color_manual(values = colors_samples)
  
pdf("../../figures/figure3.pdf", width = 12, height = 16, family = "sans")
(pA + pB) / pC / pD + 
  plot_annotation(tag_levels = "A")
dev.off()


# Running WebGestaltR takes an annoyingly long time, so I'm going to save the results
# to a file and read it back in during the plotting/formatting process. 

#WebGestaltR expects a data frame with two columns, gene name and fold change
res <- results(dds, alpha = 0.5)
res$gene <- rownames(res); rownames(res) <- NULL
res <- subset(res, select=c("gene","log2FoldChange"))
res <- as.data.frame(res)
nrow(res)

write.table(res, file = "cell_types_all_gsea.rnk", sep = "\t", quote = F, row.names = F, col.names = F)

#C8 <- suppressWarnings(WebGestaltR(enrichMethod = "GSEA",
#                                   enrichDatabaseFile = "../bulk_de/GSEA_custom_sets/c8.all.v7.5.1.symbols.gmt",
#                                   enrichDatabaseType = "genesymbol",
#                                   interestGene = res,
#                                   interestGeneType = "genesymbol",
#                                   isOutput = FALSE))
#saveRDS(C8, file = "cell_types_GSEA_results.rds")

#C8 <- readRDS("cell_types_GSEA_results.rds")

#C8 <- C8[order(C8$normalizedEnrichmentScore, decreasing = TRUE),]
#head(subset(C8, select=c("geneSet","normalizedEnrichmentScore","pValue","FDR","size")), n=10)
#tail(subset(C8, select=c("geneSet","normalizedEnrichmentScore","pValue","FDR","size")), n=10)

