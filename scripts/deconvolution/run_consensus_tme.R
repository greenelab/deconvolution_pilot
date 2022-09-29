# Consensus_TME (10.1158/0008-5472.CAN-18-3560)
# https://github.com/cansysbio/ConsensusTME

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(immunedeconv)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
source("../../config.R")

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

# Consensus_tme asks for a list of TCGA cancer types that each sample belongs to
indications <- rep("ov", length(samples))

# Run consensus_tme
res <- deconvolute_consensus_tme(bulk_matrix, indications)
res <- cbind(rownames(res), res)
colnames(res) <- c("cell_type", samples)

# Save results to text file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "consensus_tme_results.tsv", sep = "/")
write.table(res, file = text_file, sep = "\t", row.names = F, quote = F)