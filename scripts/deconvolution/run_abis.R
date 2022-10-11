# ABIS (10.1016/j.celrep.2019.01.041)
# https://github.com/giannimonaco/ABIS


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(immunedeconv)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load transcript-to-gene mapping for TPM data
tx2_gene <- fread(paste(data_path, "index/tx2gene.tsv", sep = "/"), header = F)
setnames(tx2_gene, c("Transcript", "Gene"))
gene2_symbol <- fread(paste(data_path, "index/gene2symbol.tsv", sep = "/"), header=F)
setnames(gene2_symbol, c("Gene", "Symbol"))

# Load salmon results with TPM values
bulk_matrix <- matrix()
for (i in 1:length(samples)) {
  sample_id <- samples[i]
  bulk_tmp_file <- paste(data_path, "bulk_tumors", sample_id, bulk_type,
                         "salmon/quant.sf", sep = "/")
  bulk_tmp <- fread(bulk_tmp_file)
  
  # Convert to genewise
  setnames(bulk_tmp, "Name", "Transcript")
  bulk_tmp <- right_join(tx2_gene, bulk_tmp)
  bulk_tmp <- bulk_tmp %>% group_by(Gene) %>%
    summarize(newTPM=sum(TPM))
  
  # Switch from Ensembl ID to gene names
  bulk_tmp$Gene <- gsub("\\..*", "", bulk_tmp$Gene)
  bulk_tmp <- inner_join(gene2_symbol, bulk_tmp)
  
  if (i == 1) {
    bulk_matrix <- as.matrix(bulk_tmp$newTPM)
    rownames(bulk_matrix) <- bulk_tmp$Symbol
  } else {
    bulk_matrix <- cbind(bulk_matrix, as.matrix(bulk_tmp$newTPM))
  }
}
colnames(bulk_matrix) <- samples

# Run abis
res <- deconvolute(bulk_matrix, "abis")

# The abis authors say that small negative values (>-5%) are
# expected as a result of noise and can be set to 0
res <- as.data.frame(res)
if(min(res[,-1]) < -5){
  stop("Error, found a cell type with proportion of -5% or less.
       This is a sign of strong technical or biological variability
       and this cell type should be excluded from analysis.")
}
res[res < 0] <- 0

# Switch to percentage values to allow direct comparison to other proportions
res <- cbind(cell_type=res$cell_type, res[,-1]/100)
  
# Consolidate labels based on single-cell categories
labels <- fread(paste(local_data_path,"deconvolution_input",
                      "simplified_labels.tsv", sep = "/"))
setnames(labels, "Original", "cell_type")
res <- inner_join(labels, res)

res <- res %>% group_by(Simplified) %>%
    summarize(`2251`=sum(`2251`),
              `2267`=sum(`2267`),
              `2283`=sum(`2283`),
              `2293`=sum(`2293`),
              `2380`=sum(`2380`),
              `2428`=sum(`2428`),
              `2467`=sum(`2467`),
              `2497`=sum(`2497`))
setnames(res, "Simplified", "cell_type")

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "abis_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)
