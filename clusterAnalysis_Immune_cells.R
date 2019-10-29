
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

#Plotting the immune cell propoption differences from the complete dataset

library(ggplot2)

cdf=read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/cibersort_completeteChallengeGroupDecon.csv",header = T, sep = ',')
summary(cdf)

ggplot(cdf, aes(Group, Macrophages.M1)) +
  geom_boxplot(colour = c("blue","black","green","orange","pink")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme()

#Plotting the specific challenge groups: Baseline
baseline <- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/Baseline_group.csv",header = T,sep = ',')

ggplot(baseline, aes(sample_group, Neutrophils)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

ggplot(baseline, aes(sample_group, Macrophages.M0)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

ggplot(baseline, aes(sample_group, Macrophages.M2)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

ggplot(baseline, aes(sample_group, T.cells.CD4.memory.resting)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

ggplot(baseline, aes(sample_group, Eosinophils)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

ggplot(baseline, aes(sample_group, Monocytes)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=sample_group)) + theme_classic()

# BCG challenge group
bcg <- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/BCG_challenge.csv",header = T, sep = ',')

par(mfrow=c(2,2))
ggplot(bcg, aes(Group, Monocytes)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(bcg, aes(Group, Macrophages.M0)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(bcg, aes(Group, Macrophages.M1)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(bcg, aes(Group, Macrophages.M2)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(bcg, aes(Group, T.cells.CD4.memory.resting)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(bcg, aes(Group, Eosinophils)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

#PPD challenge group

ppd <- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/PPD_challenge.csv",header = T, sep = ',')

ggplot(ppd, aes(Group, Neutrophils)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(ppd, aes(Group, Macrophages.M0)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(ppd, aes(Group, Macrophages.M1)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(ppd, aes(Group, Macrophages.M2)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(ppd, aes(Group, T.cells.CD4.memory.resting)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()

ggplot(ppd, aes(Group, Eosinophils)) +
  geom_boxplot(colour = c("blue","black","green","red")) +
  geom_jitter() +  geom_point(aes(colour=Group)) + theme_classic()


#Saline challenge group

saline <- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/nexflow_count_matrix/Saline_challenge.csv",header = T, sep = ',')

ggplot(saline, aes(Challenge_group, Neutrophils)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()

ggplot(saline, aes(Challenge_group, Macrophages.M0)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()

ggplot(saline, aes(Challenge_group, T.cells.CD4.memory.resting)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()

ggplot(saline, aes(Challenge_group, Eosinophils)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()

ggplot(saline, aes(Challenge_group, Macrophages.M2)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()

ggplot(saline, aes(Challenge_group, Monocytes)) +
  geom_boxplot(colour = c("blue","black","green","orange")) +
  geom_jitter() +  geom_point(aes(colour=Challenge_group)) + theme_classic()











































