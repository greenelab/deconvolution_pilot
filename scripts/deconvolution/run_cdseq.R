# CDSeq (10.1371/journal.pcbi.1007510)
# https://github.com/kkang7/CDSeq_R_Package

suppressPackageStartupMessages({
  library(data.table)
  library(CDSeq)
  library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples
  
# Generate bulk matrix
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

# Run CDSeq on various numbers of cell types. Note that this is a
# full deconvolution method, so no single cell reference profiles
# or marker genes are necessary
result <- CDSeq(bulk_data = bulk_matrix,
                cell_type_number = 6:15,
                mcmc_iterations = 1000)

# Save CDSeq object for later perusal
object_file <- paste(local_data_path, "deconvolution_output",
                     bulk_type, "cdseq_results_full.rds", sep = "/")
saveRDS(result, file = object_file)

# Save text version of proportion estimates, with the
# inferred number of cell types
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "cdseq_results.tsv", sep = "/")
write.table(result$estProp, file = text_file, sep = "\t")