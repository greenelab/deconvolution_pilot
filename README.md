# Effect of experimental choices on deconvolution of cancer data

## Overview

Single-cell gene expression profiling provides unique opportunities to understand tumor heterogeneity and the tumor microenvironment, but bulk profiling of tumors is still the primary population-scale analytical strategy due to cost and feasibility. Because of this, many algorithms have been developed to use single-cell profiles to deconvolve these tumors and infer their composition. However, there are a number of experimental choices made in both bulk and single-cell profiling that can bias the results of deconvolution.

In this analysis, we set out to answer the following questions:
1. Is it possible to do deconvolution with data from pooled single-cell samples as a reference, or is it necessary to spend the money to sequence all samples separately?
2. If pooling, is it better to use an antibody barcode to identify donors or to do genetic demultiplexing based on genotypes?
3. How does the dissociation process affect the composition of cell types? What biases are introduced in this step?
4. There are two common methods for increasing the amount of mRNA sequenced vs other RNA types: poly-A capture and ribosomal RNA depletion. Is one method better for samples being used for deconvolution?
5. Which existing deconvolution package is most accurate/informative for cancer data?

## Data

We've sequenced 8 high-grade serous ovarian tumors from the Penn Ovarian Cancer Research Center. For each of these samples, we have sequenced them in the following ways:
- scRNA-seq Individual (scRNA-seq, poly-A capture)
- scRNA-seq Pooled (multiplexed scRNA-seq in 2 pools of 4 samples each, poly-A capture)
- rRNA- Chunk (Bulk RNA-seq on undigested tumor chunks, rRNA depletion)
- rRNA- Dissociated (Bulk RNA-seq on dissociated cells, rRNA depletion)
- polyA+ Dissociated (Bulk RNA-seq on dissociated cells, poly-A capture)

![figure1](https://user-images.githubusercontent.com/14189222/204383085-a68437e1-33b6-42e2-a4e6-09da99040a7d.png)

The raw data (FASTQ files) is available at dbGaP: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v2.p2 (note: the dbGaP submission process is still underway. In the meantime refer to https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v1.p1 for metadata formats.)

The processed data (gene count matrices) is available at GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517

## Environment

All processing and analysis scripts were performed using the conda environment specified in environment.yml. To build and activate this environment run:

```
conda env create --file environment.yml
conda activate cancer_deconvolution
```

Once inside the conda environment, you'll need to install cellranger from the 10X website: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

To run CIBERSORTx, you'll need to obtain a license here: https://cibersortx.stanford.edu/download.php

## Adding a new method

Refer to [vignette_for_new_methods.md](https://github.com/greenelab/deconvolution_pilot/blob/main/vignette_for_new_methods.md) for step-by-step information on running a new deconvolution method through our pipeline and doing a comparative analysis.

## Citation

If you use this content, please cite:

```
Ariel A. Hippen, Dalia K. Omran, Lukas M. Weber, Euihye Jung, Ronny Drapkin, Jennifer A. Doherty, Stephanie C. Hicks, and Casey S. Greene. Performance of
computational algorithms to deconvolve heterogeneous bulk tumor tissue depends on experimental factors. bioRxiv, December 2022. 
```

View our preprint here: https://www.biorxiv.org/content/10.1101/2022.12.04.519045v1
