# deconvolution_pilot

This repository contains all the analyses for our deconvolution pilot study, which we hope will end with a set of experimental design recommendations for optimal deconvolution of cancer samples using scRNA-seq.

Questions we hope to answer:
1. Is it possible to do deconvolution with data from pooled single-cell samples as a reference, or is it necessary to spend the money to sequence all samples separately?
2. If pooling, is it better to use an antibody barcode to identify donors or to do genetic demultiplexing based on genotypes?
3. How does the dissociation process affect the composition of cell types? What biases are introduced in this step?
4. There are two common methods for increasing the amount of mRNA sequenced vs other RNA types: 3' capture (poly-A selection) and ribosomal RNA depletion. Is one method better for samples being used for deconvolution?
5. Which existing deconvolution package is most accurate/informative for cancer data?

## Data

We've sequenced 8 high-grade serous ovarian tumors from the Penn Ovarian Cancer Research Center. For each of these samples, we have sequenced them several different ways:
- scRNA-seq, 3' capture
- pooled scRNA-seq (2 pools, 4 samples each), 3' capture
- Bulk RNA-seq on tumor chunks, rRNA depletion
- Bulk RNA-seq on dissociated cells, 3' capture
- Bulk RNA-seq on dissociated cells, rRNA depletion
