# Testing additional deconvolution methods on our HGSOC data

Our deconvolution testing workflow is designed to be modular so that future users can add in other deconvolution methods and evaluate their accuracy and robustness on our data.
Follow these steps to add a new method in to our pipeline:

## Download data and index files

### Only raw read counts
If your method takes raw read counts as input, you can download the gene count matrix files from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517.

### With normalized reads
If your method requires length-normalized counts (ie transcripts per million) or if you want to compare all methods, you'll need to request to download the FastQ files from dbGaP: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v2.p2.

You will also need to download the reference files needed to run salmon. We used Gencode version 32 as a reference transcriptome (download [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz)) and the referene genome from 10X (download [here](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz)).

## Edit Snakefile and config parameters

The pipeline will be run using the snakemake workflow manager. Each step of the workflow is added as a rule (similar to a function) with explicit inputs and outputs. Which rules and in what order are calculated at runtime.
More information on how to run a snakemake pipeline is [here](https://snakemake.readthedocs.io/en/stable/).

1. Open `scripts/deconvolution/Snakefile` in your editor of choice. Look at the section labeled `# Variables`.

2. The variable called `METHODS` should have a list of deconvolution packages. You will need to modify this list in one of three ways:
  - If you only want to run your method, put `METHODS: ['yourmethodname']`
  - If you want to run your method and compare it to all raw read count methods (e.g. the ones that don't require a dbGaP download or a license from CIBERSORTx), put `METHODS: ['music', 'bisque', 'bayesprism', 'yourmethodname']`
  - If you want to run your method and compare it to all other methods, put `METHODS: ['music', 'bisque', 'bayesprism', 'epic', 'quantiseq', 'mcpcounter', 'xcell', 'consensus_tme', 'abis', 'timer', 'cibersortx', 'immucellai', 'yourmethodname']`

3. Modify the variable called `THREADS` with the desired number of cores to run on. Snakemake will use all of those threads on computationally intensive deconvolution methods, and will parallelize the less computationally intensive ones.

4. Modify all paths at the top of the Snakefile:
  - `dir_results` should be the absolute path to the directory where you want deconvolution estimates written.
  - `dir_inputs` should be the absolute path to the directory where you want intermediate files stored.
  - `dir_bulk` should be the absolute path to where bulk data was downloaded from GEO and/or dbGaP. Note that the code will expect a directory structure underneath with samples and bulk data types, e.g. `<dir_bulk>/2251/chunk_ribo/GSM6720951_bulk_dissociated_polyA_2251_STAR.tsv`
  - `dir_single` should be the absolute path to where single-cell data was downloaded from GEO. Note that the code will expect a directory structure underneath with cellranger outputs, e.g. `<dir_single>/2251/Cellranger/outs/filtered_feature_bc_matrix/matrix.mtx.gz`
  - `dir_index` should be the absolute path to any reference data downloaded.

5. Modify the paths in `config.yml` to match the paths in step 4.
  - `data_path` should be the parent directory of `dir_bulk` and `dir_single` from step 4.
  - `local_data_path` should be the parent directory of `dir_results` and `dir_inputs` from step 4.
  - `plot_path` and `figure_path` should be the absolute paths to the desired location of output plots and combined output figues. No files will be written here by snakemake, but we recommend changing them for future steps.

## Make new deconvolution script

Snakemake can be a one line bash command, but for deconvolution we recommend writing a short python or R script and using snakemake to call it.

Key components for integration with the snakemake pipeline:
- The data files can be loaded in however you see fit, but we recommend using the files made by the snakemake helper functions and stored in `dir_inputs` for compatability with other methods. 
- Be sure to specify the type of bulk data being deconvolved. In R, this can be passed directly from snakemake with the line `bulk_type <- snakemake@wildcards[['bulk_type']]`.
- The script can write any number of output files, but one of them should be titled `yourmethodname_results.tsv` with the same spelling and capitalization as in the snakemake METHODS parameter. This should be a sample by cell type matrix with the estimated cell type proportions or scores.

Below is a bare-bones example of the requirements for the deconvolution script. More detailed examples can be found in the `scripts/deconvolution` directory with the prefix `run_*`.

```
# Load required packages
library(data.table)
library(SingleCellExperiment)
library(yourMethod)
library(yaml)

# Set paramters
bulk_type <- snakemake@wildcards[['bulk_type']]
params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

# Load bulk data
bulkfile <- paste(local_data_path, "/deconvolution_input/normalized_data_", bulk_type, ".tsv", sep = "")
bulk_matrix <- fread(bulkfile)

# Load single cell data
scefile <- paste(local_data_path, "deconvolution_input", "labeled_single_cell_profile.rds", sep = "/")
sce <- readRDS(scefile)
singlecell_matrix <- assay(sce)

# Run deconvolution
deconv <- run_deconvolution_method(bulk_matrix, singlecell_matrix)

# Write to file
text_file <- paste(local_data_path, "deconvolution_output",
                   bulk_type, "yourmethod_results.tsv", sep = "/")
write.table(deconv, file=text_file)
```

## Make new rule in Snakefile

Now we will add a new rule to the Snakefile to call the new deconvolution script. The fundamental components of a snakemake rule are "input", "output", and "script". All files passed to "input" will be required to exist before the rule can be run; snakemake will look for other rules that have these files as their "output" and run them first. If the file doesn't exist and no rules generate them, snakemake will throw an error. Files in "output" will be passed to downstream rules, and eventually to the master rule `all`. Given the structure of our snakemake file, the output needs to be `yourmethodname_results.tsv`. "script" needs only to pass the name of the script written in the last step.

If your method is computationally intensive (or if you're not sure), we recommend adding in the threads parameter to give it the maximum amount of compute. 

Example rule:
```
run yourmethodname:
  input:
    "bulk_data_{bulk_type}.tsv",
    "labeled_single_cell_profile.rds"
  output:
    "yourmethodname_results.tsv"
  threads:
    THREADS
  script:
    "run_yourmethodname.R"
```

## Run snakemake
With the snakefile modified, you're now ready to run snakemake.

To check that the run looks correct (every deconvolution method is run the expected number of times, etc.) and is error-free, try:

```
snakemake --dryrun
```

It will generate the directed acyclic graph of all jobs to run and return a report of each command it will run, including all wildcard inputs.

When you're ready to run snakemake, use:

```
snakemake -j <NUMBER OF THREADS>
```

Note that running all of the existing methods on bulk and pseudo-bulk data will take several hours.

## Evaluate deconvolution output

You can visualize the robustness and accuracy of the new deconvolution results using the files in `scripts/evaluation`. All functions for loading in data and results are stored in `scripts/evaluation/evaluation_functions.R`; note that these functions may need to be lightly modified to accomodate a new user's file system.

`scripts/evaluation/get_pseudobulk_accuracy.R` returns the difference between deconvolution estimates for the pseudo-bulk sample and the true pseudo-bulk proportions, stratified by method and cell type: 
![pseudobulk_accuracy_by_cell_type](https://user-images.githubusercontent.com/14189222/204381205-c8846d7d-a881-4826-9f72-454ee59c9584.png)

`scripts/evaluation/get_accuracy.R` returns the difference between deconvolution estimates of the real bulk samples and the corresponding single cell proportions, stratified by method and cell type:
![accuracy_by_cell_type](https://user-images.githubusercontent.com/14189222/204381122-3cd9f044-12d2-4f89-a423-cfb79a1f5296.png) 

`scripts/evaluation/get_variance.R` returns the variance between deconvolution estimates across the different bulk data types in the same sample:
![variance_by_method](https://user-images.githubusercontent.com/14189222/204381298-48b7ca68-602f-4bce-bf6b-d9e6824879d3.png)

`scripts/evaluation/robustness_vs_accuracy.R` returns a two-dimensional plot of each methods' robustness and accuracy results for easy comparison: 

![accuracy_vs_robustness_real_correlation](https://user-images.githubusercontent.com/14189222/204381624-fe76b09a-3dc5-44b5-9fde-f77bf338c5fb.png)

Most of our evaluations are designed around methods that return cell type proportions for each tissue sample (e.g. "Tissue A is made up of 15% fibroblasts"). If your method returns a sample-by-sample cell type score, accuracy can be evaluated using `scripts/evaluation/get_score_correlation.R`. Be sure to add in the name of your method to the list in the `scores` variable.
