#This script is used in the visualization of the ddeconvoluted immune cells from the RNAseq data

#Data subsetting link https://dzone.com/articles/learn-r-how-extract-rows
library(magrittr)
library(tidyverse)
#Set the working directory
setwd("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization")

#Load the cell fraction data; LTBI cell fraction visualization
cell_fraction = read.table("C:/Users/Javan_Okendo/Desktop/cybersort/quantDecon/patient_cohorts/hampbel_results/hampbel_ltbi_results/quanTIseqTest_cell_fractions.txt",header = T,sep = '\t')

summary(cell_fraction)

#Boxplot for cell population comparison; this looks at the latent TB infections immune cell variations

boxplot(cell_fraction$B.cells, cell_fraction$Macrophages.M1, cell_fraction$Macrophages.M2,
        cell_fraction$Monocytes, cell_fraction$Neutrophils,cell_fraction$NK.cells,cell_fraction$T.cells.CD4,
        cell_fraction$T.cells.CD8,cell_fraction$Tregs,cell_fraction$Dendritic.cells,
        main = "Latent TB Infection: PPD,BCG,Baseline and Saline lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

#Boxplot of PPD lung challenge in the; first extract the rows containing this information

ppd_data = cell_fraction[ c(3,7,11,14), ] #extract PPD

write.table(ppd_data,"LTBI_PPD_data.txt",sep = '\t')

boxplot(ppd_data$B.cells, ppd_data$Macrophages.M1, ppd_data$Macrophages.M2,
        ppd_data$Monocytes, ppd_data$Neutrophils,ppd_data$NK.cells,ppd_data$T.cells.CD4,
        ppd_data$T.cells.CD8,ppd_data$Tregs,ppd_data$Dendritic.cells,
        main = "Latent TB Infection: PPD lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)
  
bcg_data= cell_fraction[ c(1,5,9,12), ] #Extract BCG 

write.table(bcg_data,"LTBI_BCG_data.txt",sep = '\t')

boxplot(bcg_data$B.cells, bcg_data$Macrophages.M1, bcg_data$Macrophages.M2,
        bcg_data$Monocytes, bcg_data$Neutrophils,bcg_data$NK.cells,bcg_data$T.cells.CD4,
        bcg_data$T.cells.CD8,bcg_data$Tregs,bcg_data$Dendritic.cells,
        main = "Latent TB Infection: BCG lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

baseline_data = cell_fraction[ c(2,6,10,13), ] # extract baseline data

boxplot(baseline_data$B.cells, baseline_data$Macrophages.M1, baseline_data$Macrophages.M2,
        baseline_data$Monocytes, baseline_data$Neutrophils,baseline_data$NK.cells,baseline_data$T.cells.CD4,
        baseline_data$T.cells.CD8,baseline_data$Tregs,baseline_data$Dendritic.cells,
        main = "Latent TB Infection: Baseline Sample",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

saline= cell_fraction[ c(4,8,15), ] #Extract saline data 

write.table(saline,"LTBI_saline_data.txt",sep = '\t')

boxplot(saline$B.cells, saline$Macrophages.M1, saline$Macrophages.M2,
        saline$Monocytes, saline$Neutrophils,saline$NK.cells,saline$T.cells.CD4,
        saline$T.cells.CD8,saline$Tregs,saline$Dendritic.cells,
        main = "Latent TB Infection: Saline lung challenge",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "blue",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

baseline_data_LTBI = cell_fraction[ c(2,6,10,13), ] #Extract Baseline LTBI data 

write.table(baseline_data_LTBI,"LTBI_Baseline_data.txt",sep = '\t')
#======================================================================================================================
#========================recurrent TB infection=======================================================================

#Load the data
recTB<- read.table("C:/Users/Javan_Okendo/Desktop/cybersort/quantDecon/patient_cohorts/hampbel_results/hampbel_rectb_results/quanTIseqTest_cell_fractions.txt",header=T,sep='\t')


#Extract the different patient groups 
rectb_ppd_data= recTB[ c(4,7,12), ] #rectb PPD data

write.table(rectb_ppd_data,file = "rectt_PPD_data.txt",sep = '\t',append = F,row.names = T)

#Boxplot for RECTBT 
boxplot(rectb_ppd_data$B.cells, rectb_ppd_data$Macrophages.M1, rectb_ppd_data$Macrophages.M2,
        rectb_ppd_data$Monocytes, rectb_ppd_data$Neutrophils,rectb_ppd_data$NK.cells,rectb_ppd_data$T.cells.CD4,
        rectb_ppd_data$T.cells.CD8,rectb_ppd_data$Tregs,rectb_ppd_data$Dendritic.cells,
        main = "Reccurent TB Infection: PPD lung challenge",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "blue",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

rectb_BL_data= recTB[ c(1,2,3,8,9), ] #rectb BCG data

#Boxplot for RECTBT 
boxplot(rectb_BL_data$B.cells, rectb_BL_data$Macrophages.M1, rectb_BL_data$Macrophages.M2,
        rectb_BL_data$Monocytes, rectb_BL_data$Neutrophils,rectb_BL_data$NK.cells,rectb_BL_data$T.cells.CD4,
        rectb_BL_data$T.cells.CD8,rectb_BL_data$Tregs,rectb_BL_data$Dendritic.cells,
        main = "Reccurent TB Infection: Baseline samples",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "blue",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

rectb_BCG_data= recTB[ c(6,10), ] #BCG

#Boxplot for RECTBT 
boxplot(rectb_BCG_data$B.cells, rectb_BCG_data$Macrophages.M1, rectb_BCG_data$Macrophages.M2,
        rectb_BCG_data$Monocytes, rectb_BCG_data$Neutrophils,rectb_BCG_data$NK.cells,rectb_BCG_data$T.cells.CD4,
        rectb_BCG_data$T.cells.CD8,rectb_BCG_data$Tregs,rectb_BCG_data$Dendritic.cells,
        main = "Reccurent TB Infection: BCG lung challenge immune profile ",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "blue",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)

rectb_sal_data= recTB[ c(5,13), ] #sal

#Boxplot for RECTBT 
boxplot(rectb_sal_data$B.cells, rectb_sal_data$Macrophages.M1, rectb_sal_data$Macrophages.M2,
        rectb_sal_data$Monocytes, rectb_sal_data$Neutrophils,rectb_sal_data$NK.cells,rectb_sal_data$T.cells.CD4,
        rectb_sal_data$T.cells.CD8,rectb_sal_data$Tregs,rectb_sal_data$Dendritic.cells,
        main = "Reccurent TB Infection: Saline lung challenge immune profile ",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "blue",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)


#=========================================================================================================================
#==ABSI cell deconvolution from the human lung samples following different challenges with different molecules

Abis_immune_cells <- read.table("C:/Users/Javan_Okendo/Desktop/cybersort/quantDecon/patient_cohorts/hampbel_results/hampbel_rectb_results/quanTIseqTest_cell_fractions.txt",sep = '\t',header = T)

#Principal component analysis
library("factoextra")
# before the PCA analysis
res.pca <- prcomp(Abis_immune_cells[, -1],  scale = TRUE)

# Default plot
fviz_pca_ind(res.pca)


# Use points only
fviz_pca_ind(res.pca, geom="point",pointsize = 4)

# Change the theme and use only points
fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)+ theme_minimal()

# Color individuals by groups
p <- fviz_pca_ind(res.pca, label="none",addEllipses=TRUE, 
                  ellipse.level=0.95,habillage=Abis_immune_cells$Sample, 
                  title="Reccurent TB patient group immune profile after PPD,BCG,Saline and Baseline lung challenge")

print(p)

#==================================================================================
#=======Previous Tuberculosis======================================================

prevTB_cohort=read.table("C:/Users/Javan_Okendo/Desktop/cybersort/quantDecon/patient_cohorts/prevTB_decon_results/quanTIseqTest_cell_fractions.txt",header=T,sep='\t')

View(prevTB_cohort)


boxplot(prevTB_cohort$B.cells, prevTB_cohort$Macrophages.M1, prevTB_cohort$Macrophages.M2,
        prevTB_cohort$Monocytes, prevTB_cohort$Neutrophils,prevTB_cohort$NK.cells,prevTB_cohort$T.cells.CD4,
        prevTB_cohort$T.cells.CD8,prevTB_cohort$Tregs,prevTB_cohort$Dendritic.cells,
        main = "Previous TB Infection: PPD,BCG,Baseline and Saline lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white")
)
#=============================================================================================
#Extract the patient group according to the lung challenge
prevtb_BCG_data= prevTB_cohort[ c(1,5,10,13,17,21,25,29,34,38,42,50,56,60), ] #BCG challenge

#Boxplot
boxplot(prevtb_BCG_data$B.cells, prevtb_BCG_data$Macrophages.M1, prevtb_BCG_data$Macrophages.M2,
        prevtb_BCG_data$Monocytes, prevtb_BCG_data$Neutrophils,prevtb_BCG_data$NK.cells,prevtb_BCG_data$T.cells.CD4,
        prevtb_BCG_data$T.cells.CD8,prevtb_BCG_data$Tregs,prevtb_BCG_data$Dendritic.cells,
        main = "Previous TB Infection: BCG lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white"))


prevtb_PPD_data= prevTB_cohort[ c(3,7,12,15,19,23,27,31,36,40,44,48,52,58,62), ] #PPD challenge

#Boxplot
boxplot(prevtb_PPD_data$B.cells, prevtb_PPD_data$Macrophages.M1, prevtb_PPD_data$Macrophages.M2,
        prevtb_PPD_data$Monocytes, prevtb_PPD_data$Neutrophils,prevtb_PPD_data$NK.cells,prevtb_PPD_data$T.cells.CD4,
        prevtb_PPD_data$T.cells.CD8,prevtb_PPD_data$Tregs,prevtb_PPD_data$Dendritic.cells,
        main = "Previous TB Infection: PPD lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white"))

prevtb_BL_data= prevTB_cohort[ c(2,6,9,11,14,18,22,26,30,33,35,39,43,46,47,51,54,55,57,61), ] #Baseline challenge

#Boxplot
boxplot(prevtb_BL_data$B.cells, prevtb_BL_data$Macrophages.M1, prevtb_BL_data$Macrophages.M2,
        prevtb_BL_data$Monocytes, prevtb_BL_data$Neutrophils,prevtb_BL_data$NK.cells,prevtb_BL_data$T.cells.CD4,
        prevtb_BL_data$T.cells.CD8,prevtb_BL_data$Tregs,prevtb_BL_data$Dendritic.cells,
        main = "Previous TB Infection: Baseline samples",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white"))

prevtb_sal_data= prevTB_cohort[ c(4,8,16,20,24,28,32,37,41,45,49,53,59,63), ] #Saline challeng

#Boxplot
boxplot(prevtb_sal_data$B.cells, prevtb_sal_data$Macrophages.M1, prevtb_sal_data$Macrophages.M2,
        prevtb_sal_data$Monocytes, prevtb_sal_data$Neutrophils,prevtb_sal_data$NK.cells,prevtb_sal_data$T.cells.CD4,
        prevtb_sal_data$T.cells.CD8,prevtb_sal_data$Tregs,prevtb_sal_data$Dendritic.cells,
        main = "Previous TB Infection: Saline lung challenge immune profile",
        names = c("B_cells","Macrophages.M1","Macrophages.M2","Monocytes","Neutrophils","NK_cells","T_cells_CD4","T_cells_CD8",
                  "Tregs","Dendritic cells"),
        las = 1.5,
        border = "brown",
        horizontal = F,
        notch = F,
        xlab="Immune Cell type",
        ylab="Quantified immune Expression level",
        col = ("white"))


#Principal component analysis
library("factoextra")
# before the PCA analysis
res.pca <- prcomp(prevTB_cohort[, -1],  scale = TRUE)

# Default plot
fviz_pca_ind(res.pca)


# Use points only
fviz_pca_ind(res.pca, geom="point",pointsize = 4,col.ind = "contrib",
             title="Previous TB patient group immune profile after PPD,BCG,Saline and Baseline lung challenge")

# Change the theme and use only points
fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") +
  scale_color_gradient2(low="white", mid="blue",
  high="red", midpoint=0.6)+ theme_minimal()

# Color individuals by groups
p <- fviz_pca_ind(res.pca, label="none",addEllipses=F, 
                  ellipse.level=0.95,habillage=prevTB_cohort$Sample, 
                  title="Reccurent TB patient group immune profile after PPD,BCG,Saline and Baseline lung challenge")

print(p)

#==================easyboxplot=======================
setwd("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization/complete_Tx_data_decon")

library(devtools)
library(easyGgplot2)
#install_github("kassambara/easyGgplot2", dependencies = T)

df=read.table("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization/complete_Tx_data_decon/quanTIseqTest_cell_fractions.txt",header=T,sep='\t')

head(df, 6)


library(Rtsne)
## Learning t-SNE Plotting
## Load dataset
IR <- df # Loading the immune  dataset into a object called IR

## Split IR into two objects: 1) containing measurements 2) containing species type
IR_data <- IR[ ,2:11] # We are sub-setting IR object such as to include 'all rows' and columns 3 to 12.
IR_species <- IR[ ,2] # We are sub-setting IR object such as to include 'all rows' and column 2.

## Load the t-SNE library
library(Rtsne)

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(IR_data, perplexity=5, check_duplicates = FALSE,verbose=T,max_iter=500) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "blue", pch = 19, cex = 1.5,xlab = "t-SNE 1",ylab = "t-SNE 2",
     main = "complete Tx data: PPD,BCG, Baseline and Saline Lung challenge group") # Plotting the first image
plot(tsne_results$Y, col = "black", bg= IR_species, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "complete Tx data: PPD,BCG, Baseline and Saline Lung challenge") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=IR$Group)



#============barplot

#. 1 LTBI
ltbi=read.table("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization/complete_Tx_data_decon/LTBI_combined_cellfraction.txt",header = T,sep = '\t')

head(ltbi,4)

# Create the data for the chart
H <- ltbi[,2:16]

h = as.matrix(H)

M = c("B.cells","M1","M2","Monocytes","Neutrophils","NK.cells","CD4","CD8","Tregs","DC")

k=c("T004_BCG","T006_BCG","T016_BCG","T042_BCG",
    "T004_PPD","T006_PPD","T016_PPD","T042_PPD",
    "T004_BL","T006_BL","T016_BL","T042_BL","T004_sal",
    "T006_sal","T042_sal")


color= c("green","orange","brown","red","blue","yellow","Lavender","Black","Cyan","Chocolate","Gold")
# Plot the bar chart 
barplot(h,names.arg=k,xlab="Immune cell type",ylab="Expression proportion",col= color,
        main="Propoprtion of Immune cells per sample in LTBI patient Group ",srt=45,border="red",las=1.5)

# Add the legend to the chart
legend("topleft", M, cex = 1.3, fill = color)
# Save the file
#======================================================

#2 recurrent TB
rectb=read.table("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization/complete_Tx_data_decon/rectTB_transposed.txt",header = T,sep = '\t')

head(rectb,4)

# Create the data for the chart
sub_1 <- rectb[,2:14]

h = as.matrix(sub_1)

M = c("B.cells","M1","M2","Monocytes","Neutrophils","NK.cells","CD4","CD8","Tregs","DC")

k=c("T051_BCG","T054_BCG","T051_PPD","T051_PPD","T054_PPD","T049_BL","T050_BL","T051_BL","T052_BL","T053_BL","T054_BL","T051_sal","T054_sal")


color= c("green","orange","brown","red","blue","yellow","Lavender","Black","Cyan","Chocolate")
# Plot the bar chart 
barplot(h,names.arg=k,xlab="Immune cell type",ylab="Expression proportion",col= color,
        main="Propoprtion of Immune cells per sample in recurrent TB patient Group ",border="red",las=1.5)

# Add the legend to the chart
legend("topleft", M, cex = 1.3, fill = color)
# Save the file

#=========================================================

#3. previuos TB

#2 recurrent TB
prevTB=read.table("C:/Users/Javan_Okendo/Desktop/cybersort/prism_visualization/complete_Tx_data_decon/prevTB_transposed.txt",header = T,sep = '\t')

head(prevTB,4)

# Create the data for the chart
sub_1 <- prevTB[,2:64]

h = as.matrix(sub_1)

M = c("B.cells","M1","M2","Monocytes","Neutrophils","NK.cells","CD4","CD8","Tregs","DC")

k=c("T007_BCG","T009_BCG","T011_BCG","T013_BCG","T014_BCG","T015_BCG","T025_BCG","T026_BCG",
    "T031_BCG","T032_BCG","T033_BCG","T039_BCG","T032_BCG","T047_BCG","T048_BCG","T007_PPD","T009_PPD",
    "T011_PPD","T013_PPD","T014_PPD","T015_PPD","T025_PPD","T026_PPD","T031_PPD","T032_PPD","T033_PPD","T038_PPD","T039_PPD","T047_PPD","T048_PPD","T007_BL",
    "T009_BL","T010_BL","T011_BL","T013_BL","T014_BL","T015_BL","T025_BL","T026_BL","T030_BL","T032_BL","T033_BL","T034_BL","T038_BL","T039_BL",
    "T045_BL","T046_BL","T047_BL","T048_BL","T007_sal","T009_sal","T013_sal","T014_sal","T015_sal","T025_sal",
    "T026_sal","T031_sal","T032_sal","T033_sal","T038_sal","T039_sal","T047_sal","T048_sal")


color= c("green","orange","brown","red","blue","yellow","Lavender","Black","Cyan","Chocolate")
# Plot the bar chart 
barplot(h,names.arg=k,xlab="Immune cell type",ylab="Expression proportion",col= color,
        main="Propoprtion of Immune cells per sample in recurrent TB patient Group ",border="red",las=1.5)

# Add the legend to the chart
legend("topleft", M, cex = 1.3, fill = color)
# Save the file
































































