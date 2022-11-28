# Testing additional deconvolution methods on our data

Our deconvolution testing workflow is designed to be modular so that future users can add in other deconvolution methods and evaluate their accuracy and robustness on our data.
Follow these steps to add a new method in to our pipeline:

## Download data and index files

### Only raw read counts
If your method takes raw read counts as input, you can download the gene count matrix files from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517 

### With normalized reads
If your method requires length-normalized counts (ie transcripts per million) or if you want to compare all methods, you'll need to request to download the FastQ files from dbGaP: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v2.p2.

You will also need to download the reference files needed to run salmon. We used Gencode version 32 as a reference transcriptome (download [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz)) and the referene genome from 10X (download [here](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz)).

## Edit Snakefile and config parameters

The pipeline will be run using the snakemake workflow manager. Each step of the workflow is added as a rule (similar to a function) with explicit inputs and outputs. Which rules and in what order are calculated at runtime.
More information on how to run a snakemake pipeline is [here](https://snakemake.readthedocs.io/en/stable/).
