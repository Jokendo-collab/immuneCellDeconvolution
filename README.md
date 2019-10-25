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
### Load the count matrix dataframe
`Countdata = read.delim("fullChallengGroups_gene_counts.txt",header = T,row.names = 1)`

### check if the row names in metadata matches with the column names in the count matrix data
`ncol(Countdata) == nrow(metadata)`

### Load the metadata file
```
metadata = read.csv("metadata.txt",header = T,sep = '\t')

rownames(metadata) <- colnames(Countdata) #match the rownames in metadata with colnames in the count matrix data

rownames(metadata) == colnames(Countdata)

ddsMat <- DESeqDataSetFromMatrix(countData = Countdata,
                                 colData = metadata,
                                 design = ~Replicate)


```
### Find differential expressed genes
```
# Run DESEq2
ddsMat <- DESeq(ddsMat)
```

_Get basic statisics about the number of significant genes_
### Get results from testing with FDR adjust pvalues
`results <- results(ddsMat, pAdjustMethod = "BH", alpha = 0.05)`

### Generate summary of testing. 
`summary(results)`

### Check directionality of fold change
`mcols(results, use.names = T)`

### - - - - - - - - - - - - - 
### Gene annotation
### - - - - - - - - - - - - - 
```
library(AnnotationDbi)
library(org.Hs.eg.db)
```
## Add gene full name
```
results$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05)
head(results_sig)

# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), 
            file = 'normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(results), 
            file = "results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

```
### - - - - - - - - - - - - - 
### PCA plot
### - - - - - - - - - - - - - 
### Convert all samples to rlog
`ddsMat_rlog <- rlong(ddsMat, blind = FALSE)`
`ddsMat_rlog <- vst(ddsMat, blind = FALSE)`
### Plot PCA by column variable
```
plotPCA(ddsMat_rlog, intgroup = "Replicate", ntop = 500) +
  theme_bw() +
  ggsave('challeng_group_plot.png')

plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) +
  theme_bw() +
  ggsave('group_plot.png')
  
# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Replicate", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 3) + # Increase point size
  scale_y_continuous(limits = c(-25, 25)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") 
```
### - - - - - - - - - - - - - 
### Heatmap plot
### - - - - - - - - - - - - - 
### Load libraries
```
library(pheatmap) 
library(RColorBrewer) 

# Convert all samples to rlog
ddsMat_rlog <- vst(ddsMat, blind = FALSE)

# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:50, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group),
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$SampleID
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c("PTB" = "blue", "RTB" = "black","LTBI" = "green","STI" = "red"),
  Replicate = c(Baseline = "red", PPD = "green", BCG= "yellow", Saline = "blue")
)
```
### Make Heatmap with pheatmap function
_See more in documentation for customization_
```
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(25), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 5, # Make the cells wider
         show_colnames = T)
```
### - - - - - - - - - - - - - 
### Volcano plot
### - - - - - - - - - - - - - 
```
# Load libraries
library(ggplot2)
library(RColorBrewer)

# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue < 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval < 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("PPD" / "BCG"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

```
### - - - - - - - - - - - - - 
### MA plot
### - - - - - - - - - - - - - 
`plotMA(results, ylim = c(-5, 5))`

#plot dispersion
plotDispEsts(ddsMat)

### - - - - - - - - - - - - - 
### Single gene plot
### - - - - - - - - - - - - - 
### Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

### Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

### Plot single gene
plotCounts(ddsMat, gene = top_gene, intgroup = "Group",col=rainbow(4))
### - - - - - - - - - - - - - - -
# Pathway analysis of DE genes
### - - - - - - - - - - - - - - -

### Load required libraries
library(clusterProfiler)
library(ReactomePA)
library(KEGGREST)
library(DOSE)
library(org.Mm.eg.db)
library(pathview)

# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)

# - - - - - - - - - - - - -
# Enrich with KEGG database
# - - - - - - - - - - - - -
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'human',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)


# - - - - - - - - - - - - -
# Enrich with ReactomeDB
# - - - - - - - - - - - - -
reactome_enrich <- enrichPathway(gene = names(gene_matrix), 
                                 organism = 'human', 
                                 readable = TRUE,
                                 pvalueCutoff = 0.05)
# Get table of results
head(summary(reactome_enrich))

# Plot results
barplot(reactome_enrich, 
        drop = TRUE, 
        showCategory=10, 
        title = "ReactomeDB Enrichment Pathways",
        font.size = 8)


# - - - - - - - - - - - - -
# Enrich with GO
# - - - - - - - - - - - - -
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological process",
        font.size = 8)

# - - - - - - - - - - - - -
# Plot specific KEGG pathways
# (with fold change) 
# - - - - - - - - - - - - -
# pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "04070", 
         species = "human", 
 ```
         map.symbol = T)
