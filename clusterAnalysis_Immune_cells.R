
#Load the dataframe
compImmuneData<- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/complete_BAL_countmatrix/complete_sample_group_immune_deconvolution.csv",header = T,sep = ',')

head(compImmuneData) #Check the first few lines of the data

BCG_challenge <- compImmuneData[grep("BCG",compImmuneData$Group), ] #Extract the BCG patient Cohort

PPD_challenge <- compImmuneData[grep("PPD",compImmuneData$Group), ]

Baseline_group <- compImmuneData[grep("Baseline",compImmuneData$Group), ]

Saline_challenge <- compImmuneData[grep("Saline",compImmuneData$Group), ]

#==============Cluster analysis of different challenge groups

library(Rtsne)
## Learning t-SNE Plotting
## Load dataset
df <- BCG_challenge # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=5000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "BCG challenge group: PrevTB,LTBI & RecTB") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Input.Sample)

#PPD challenge group t-SNE plot

## Learning t-SNE Plotting
## Load dataset
df <- PPD_challenge # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=5000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "PPD challenge group: PrevTB,LTBI & RecTB") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Input.Sample)

#Baseline patient group

## Learning t-SNE Plotting
## Load dataset
df <- Baseline_group # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=5000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "Baseline group: PrevTB,LTBI & RecTB") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Input.Sample)

#Saline challenge group
## Learning t-SNE Plotting
## Load dataset
df <- Saline_challenge # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=5000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "Saline challenge group: PrevTB,LTBI & RecTB") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Input.Sample)

#Complete patient group

## Learning t-SNE Plotting
## Load dataset
df <- compImmuneData # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=1, check_duplicates = FALSE,verbose=T,max_iter=50000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "Complete challenge group: PrevTB,LTBI & RecTB") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Group)


#PCA plot
library("FactoMineR")
library("factoextra")

res.pca <- PCA(cell_proportions, graph = FALSE)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#Plot variable type
fviz_pca_var(res.pca, col.var = "black")

fviz_pca_ind(res.pca,  label="none", habillage=compImmuneData$Group,addEllipses=F, ellipse.level=0.95)

#K-means clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#=======================================================================================
#Plotting complexHeatmap
# This is the reeanalysis script including sterilizing immunity

setwd("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/")
library(ComplexHeatmap)
##Load the data
cdf=read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/cibersort_completeteChallengeGroupDecon.csv",header = T, sep = ',')
summary(cdf)


BCG_challenge <- cdf[grep("_BCG",cdf$Input.Sample), ] #Extract the BCG patient Cohort

##write.csv(BCG_challenge,file = "BCG_challenge.csv")

PPD_challenge <- cdf[grep("_PPD",cdf$Input.Sample), ]

##write.csv(PPD_challenge,file = "PPD_challenge.csv")

Baseline_group <- cdf[grep("_BL",cdf$Input.Sample), ]

##write.csv(Baseline_group,file = "Baseline_group.csv")

Saline_challenge <- cdf[grep("_Sal",cdf$Input.Sample), ]

##write.csv(Saline_challenge,file = "Saline_challenge.csv")

sterilizing_immu <- cdf[grep("_SI_",cdf$Input.Sample), ]

##write.csv(sterilizing_immu,file = "sterilizing_immu.csv")


## Learning t-SNE Plotting
## Load dataset
library(Rtsne)
df <- cdf # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing immune cell proportions 2) containing challenge groups
cell_proportions <- cdf[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 2 to 22.
challenge_groups <- cdf[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(cell_proportions, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=1000) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
plot(tsne_results$Y, col = "black", bg= challenge_groups, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "Complete challenge group: PrevTB,LTBI, RecTB $ Sterilizing Immunity") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=cdf$Group)
#===============================================================================================

















































