# Introduction
The cell type composition of many biological tissues varies widely across samples. Such sample heterogeneity hampers efforts to probe the role of each cell type in the tissue. The use of in silico approaches to complement the flow cytometer is of great importance as it will help in `"picking"` the immune cells at the transcription level which otherwise could not be "seen" at the proteomic level. The RNAset deconvolution provides more power in making a conclusion about a given immune cell population in a specific tissue or microenvironment.
## The workflow shows the summary of the mixture matrix generation from the RNA fastq files
![w2](https://user-images.githubusercontent.com/26459707/66849989-35d86980-ef78-11e9-9971-fac9fdd9e1a9.png)
# The pipeline does the following
The workflow processes raw data from FastQ inputs (`FastQC`, `Trim Galore!`), aligns the reads using STAR sligner, generates gene counts with the featureCounts tooland performs extensive quality-control on the results (`RSeQC`, `dupRadar`, `Preseq`, `edgeR`, `MultiQC`).
### Quality control
This is used to assess the quality of the reads from the sequencing facility
### Treaming and alignment of the reads to the reference genome
To remove the adapter sequences and also do the quality control on the trimmed reads
### Read summarization using subreads
Quantify the number of genes mapped to the reference genome
### Immune cell fraction estimation using `CIBERSORT` pipeline
Estimation of the Immune cell fractions present from the RNAseq data. It uses unsupervised machine learning approach to estimate the immune cells from our bulk sample**. **Scirpts for `CIBERSORT.R` can be accessed from (https://cibersort.stanford.edu/) upon an request from `CIBERSORT` team.

# Differential Analysis of the count matrix
```
#set the working directory
setwd("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix")
list.files()
# Load the following Bioconductor libraries and if you don't have them installed then you can do so from the Bioconductor 
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(ggsci)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(KEGGREST)
library(gage)
library(sqldf)
```
