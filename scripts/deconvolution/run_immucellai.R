# ImmuCellAI (10.1002/advs.201902880)
# https://github.com/lydiaMyr/ImmuCellAI

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ImmuCellAI)
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

# group_tag is 0 because we're not trying to compare between different groups
# response_tag is 0 because we don't need to predict response to ICB therapy
res <- ImmuCellAI_new(bulk_matrix, "rnaseq", group_tag = 0, response_tag = 0)

# Save full object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "immucellai_results_full.rds", sep = "/")
saveRDS(res, file = object_file)

# Reformat output data to match style of other deconvolution methods
out <- t(res$Sample_abundance)
out <- cbind(rownames(out),out)
colnames(out) <- c("cell_type", samples)

outfile <- paste(local_data_path, "deconvolution_output", bulk_type,
                 "immucellai_results.tsv", sep = "/")
write.table(out, outfile, sep = "\t", quote=F, row.names = F)