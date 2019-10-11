# Introduction
Integration of heterogeneous and voluminous data from proteomics, transcriptomics, Immunological and clinical research constitutes not only a fundamental problem but a real hurdle in the extraction of valuable information from these omics data sets. The exponential increase of novel omics technologies such as LC-MS/MS and the generation of high-resolution data from large consortia projects generates heterogenous and big data sets. These big data promote research but at the same time, new methodologies are lacking in sophisticated tools to analyse these complex data sets. Here we present horizontal data integration approach.
## The workflow shows the summary of the mixture matrix generation from the RNA fastq files
![w2](https://user-images.githubusercontent.com/26459707/66654239-14a31080-ec3a-11e9-9e6a-83b648c7b6fd.png)
# The pipeline does the following:
### Quality control
This is used to assess the quality of the reads from the sequencing facility
### Treaming and alignment of the reads to the reference genome
To remove the adapter sequences and also do the quality control on the trimmed reads
### Read summarization using subreads
Quantify the number of genes mapped to the reference genome 
### Immune cell fraction estimation using CIBERSORT pipeline
Estimation of the Immune cell fractions present from the RNAseq data. It uses unsupervised machine learning approach to estimate the immune cells from our bulk sample

