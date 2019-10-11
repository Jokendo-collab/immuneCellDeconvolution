library(devtools)
##install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(UpSetR)

cellPop = read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/cybersort_source_code/TB_hart_deconPDF/complete_sample_group_immune_deconvolution.csv", 
                  header = TRUE, sep = ",")
head(cellPop)

Heatmap(cellPop[2:23])

nrow(cellPop)

ncol(cellPop)

#============================
#Heatmap analysis of the different cell fractions
library(ggplot2)
library(RColorBrewer)
library(heatmap.plus)
#Previous TB Patient groups

#Load the dataframe
prevTB<- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/cybersort_source_code/TB_hart_deconPDF/prevTB.csv",header = T,sep = ',')

head(prevTB) #Check the first few lines of the data

BCG <- prevTB[grep("_BCG",prevTB$Input.Sample), ] #Extract the BCG patient Cohort

PPD <- prevTB[grep("_PPD",prevTB$Input.Sample), ]

Baseline <- prevTB[grep("_BL",prevTB$Input.Sample), ]

Saline <- prevTB[grep("_sal",prevTB$Input.Sample), ]

#==========================PCA============================

df = read.csv("complete_sample_group_immune_deconvolution.csv",header = T,sep = ',')
library("FactoMineR")
library("factoextra")
df2=df[,-1:-2]
PCA(df2, scale.unit = TRUE, ncp = 5, graph = F)

#===================correlation
res.pca <- PCA(df2, graph = FALSE)
print(res.pca)

#Visualize and interprate the data
fviz_pca_biplot(res.pca)

#Screeplot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#Correlation between the immune cells in BALF
fviz_pca_var(res.pca, col.var = "red",repel = T)

#Get the PCA variables
var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

#correlation plot
library("corrplot")
corrplot(var$cos2, is.corr=F)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#Run the PCA function
iris.pca <- PCA(df2, graph = FALSE)

fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df$Group, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#000066"),
             addEllipses = F, # Concentration ellipses
             legend.title = "Challenge Groups"
)

#Naming the individual points
ind.p <- fviz_pca_ind(iris.pca, geom = "point", col.ind = df$Group)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Baseline patient challenge groups",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Species", legend.position = "top",
              ggtheme = theme_light(), palette = "jco"
)


#Add more aestethics
fviz_pca_biplot(iris.pca,repel = TRUE,
                # Individuals
                geom.ind = "point",
                fill.ind = df$Group, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = F,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Group", color = "Contrib",
                                    alpha = "Contrib")
)
#=====================End of PCA========================================================


setwd("C:/Users/Javan_Okendo/Desktop/cybersort/complete_BAL_countmatrix/")
total_samples <- read.csv("C:/Users/Javan_Okendo/Desktop/cybersort/complete_BAL_countmatrix/complete_sample_group_immune_deconvolution.csv",header = T,sep = ',')

#Extract the different patient groups for t-SNE analysis

LTBI <- total_samples[grep("_LTBI",total_samples$Input.Sample), ] #LTBI patient group
write.csv(LTBI, file = "LTBI_patient_group_subset.csv")

prevTB <- total_samples[grep("_prevTB",total_samples$Input.Sample), ] #PPD patient group
write.csv(prevTB,file = "prevTB_patient_group_subset.csv")

rectb <- total_samples[grep("_recTB_",total_samples$Input.Sample), ] #saline patient group

write.csv(rectb,file = "recTB_patient_group_subset.csv")

#==============Latent TB infection
BL_group<- read.csv("recTB_patient_group_subset.csv",sep = ',',header = T)

library(Rtsne)
## Learning t-SNE Plotting
## Load dataset
df <- BL_group # Loading the immune  dataset into a object called IR

## Split df into two objects: 1) containing measurements 2) containing species type
IR_data <- df[ ,2:22] # We are sub-setting df object such as to include 'all rows' and columns 3 to 12.
IR_species <- df[ ,2] # We are sub-setting df object such as to include 'all rows' and column 2.

## Load the t-SNE library
library(Rtsne)

## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(IR_data, perplexity=2, check_duplicates = FALSE,verbose=T,max_iter=500) # You can change the value of perplexity and see how the plot changes

## Generate the t_SNE plot
par(mfrow=c(1,1)) # To plot two images side-by-side
plot(tsne_results$Y, col = "black", bg= IR_species, pch = 21, cex = 1.5,xlab = "t-SNE 1",
     ylab = "t-SNE 2",main = "Recurrent TB Group: PPD,BCG, Baseline and Saline Lung challenge") # Second plot: Color the plot by the real species type (bg= IR_species)
text(tsne_results$Y, labels=df$Group)

#t-SNE plot with ggplot package
library(ggplot2)
tsne_plot <- data.frame(x = tsne_results$Y[,1], y = tsne_results$Y[,2], Challenge_Group = df$Group)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=Challenge_Group))+theme_bw()+
xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("Previous TB: PPD,BCG, Baseline $ Saline lung Challenge")+
theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))




























































