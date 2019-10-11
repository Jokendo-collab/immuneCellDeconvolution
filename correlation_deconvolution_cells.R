setwd("C:/Users/Javan_Okendo/Desktop/cybersort/decon_statistical_test")

library(dplyr)

prevTB <- read.csv("prevTB.csv",header = T,sep = ',')

head(prevTB)

attach(prevTB)

# Show the g levels
levels(prevTB$g)


#Previous TB
# Box plots
# Plot weight by g and color by g
library("ggpubr")
colnames(prevTB)

#B.cell naive
par(mfrow=c(4,4)) # To plot two images side-by-side
ggboxplot(prevTB, x = "Group", y = "Neutrophils", 
          color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Neutrophils", xlab = "Challenge gs",title = "Neutrophils")


ggline(prevTB, x = "Group", y = "T.cells.regulatory..Tregs.", 
       add = c("mean_se", "jitter"), 
       order = c("Baseline","BCG","PPD","Saline"),
       ylab = "T.cells.regulatory..Tregs.", xlab = "Challenge g",title = "Mast.cells.resting")


#computation of Kruskal-Wallis test
kruskal.test(Macrophages.M0 ~ Group, data = prevTB)

#Multiple pairwise-comparison between gs
pairwise.wilcox.test(prevTB$Macrophages.M0, prevTB$Group,
                     p.adjust.method = "BH")

#Recurrent TB patient g

rectb <- read.csv("reccurentTB_cibersort_results.csv",header = T, sep = ',')

attach(rectb)
colnames(rectb)

# this will be done interactively for the 22 immune cells per g
ggboxplot(rectb, x = "g", y = "Neutrophils", 
          color = "g", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Neutrophils", xlab = "Challenge gs",title = "Neutrophils")


#Latent TB infection
LTBI <- read.csv("LTBI_cibersort_results.csv",header = T,sep = ',')
attach(LTBI)

colnames(LTBI)

# this will be done interactively for the 22 immune cells per g
ggboxplot(LTBI, x = "challeng_group", y = "B.cells.naive", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "B.cells.naive", xlab = "Challenge gs",title = "B.cells.naive")
#=======================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "B.cells.memory", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "B.cells.memory", xlab = "Challenge gs",title = "B.cells.memory")
#========================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "Plasma.cells", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Plasma.cells", xlab = "Challenge gs",title = "Plasma.cells")

#==============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.CD8", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.CD8", xlab = "Challenge gs",title = "T.cells.CD8")
#===============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.CD4.naive", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.CD4.naive", xlab = "Challenge gs",title = "T.cells.CD4.naive")
# ===============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.CD4.memory.resting", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.CD4.memory.resting", xlab = "Challenge gs",title = "T.cells.CD4.memory.resting")
#================================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.CD4.memory.activated", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.CD4.memory.activated", xlab = "Challenge gs",title = "T.cells.CD4.memory.activated")

# ===============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.follicular.helper", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.follicular.helper", xlab = "Challenge gs",title = "T.cells.follicular.helper")
#================================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "T.cells.regulatory..Tregs.", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.regulatory..Tregs.", xlab = "Challenge gs",title = "T.cells.regulatory..Tregs.")
# ===============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "T.cells.gamma.delta", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "T.cells.gamma.delta", xlab = "Challenge gs",title = "T.cells.gamma.delta")
# ===============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "NK.cells.resting", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "NK.cells.resting", xlab = "Challenge gs",title = "NK.cells.resting")
# ===============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "NK.cells.activated", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "NK.cells.activated", xlab = "Challenge gs",title = "NK.cells.activated")
# ===============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Monocytes", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Monocytes", xlab = "Challenge gs",title = "Monocytes")
# ===============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "Macrophages.M0", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Macrophages.M0", xlab = "Challenge gs",title = "Macrophages.M0")
# ===============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Macrophages.M1", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Macrophages.M1", xlab = "Challenge gs",title = "Macrophages.M1")
# ==============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Macrophages.M2", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Macrophages.M2", xlab = "Challenge gs",title = "Macrophages.M2")
# ==============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Dendritic.cells.resting", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Dendritic.cells.resting", xlab = "Challenge gs",title = "Dendritic.cells.resting")
# ==============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Dendritic.cells.activated", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Dendritic.cells.activated", xlab = "Challenge gs",title = "Dendritic.cells.activated")
# ==============================================================================================================

ggboxplot(LTBI, x = "challeng_group", y = "Mast.cells.resting", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Mast.cells.resting", xlab = "Challenge gs",title = "Mast.cells.resting")
# =============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "Mast.cells.activated", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Mast.cells.activated", xlab = "Challenge gs",title = "Mast.cells.activated")

#==============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "Eosinophils", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Mast.cells.activated", xlab = "Challenge gs",title = "Eosinophils")

#==============================================================================================================
ggboxplot(LTBI, x = "challeng_group", y = "Neutrophils", 
          color = "challeng_group", palette = c("#00AFBB", "#E7B800", "#FC4E07","#202020"),
          order = c("Baseline","BCG","PPD","Saline"),
          ylab = "Neutrophils", xlab = "Challenge gs",title = "Neutrophils")

#Assessing the distribution of data
par(mfrow=c(5,5))
#Previous TB
hist(prevTB$B.cells.naive,main = "B.cells.naive")
hist(prevTB$B.cells.memory,main = "B.cells.memory")
hist(prevTB$Plasma.cells,main = "Plasma.cells")
hist(prevTB$T.cells.CD8,main = "T.cells.CD8")
hist(prevTB$T.cells.CD4.naive,main = "T.cells.CD4.naive")
hist(prevTB$T.cells.CD4.memory.resting,main = "T.cells.CD4.memory.resting")
hist(prevTB$T.cells.CD4.memory.activated,main = "T.cells.CD4.memory.activated")
hist(prevTB$T.cells.follicular.helper,main = "T.cells.follicular.helper")
hist(prevTB$T.cells.regulatory..Tregs.,main = "T.cells.regulatory..Tregs.")
hist(prevTB$T.cells.gamma.delta,main = "T.cells.gamma.delta")
hist(prevTB$NK.cells.resting,main = "NK.cells.resting")
hist(prevTB$NK.cells.activated,main = "NK.cells.activated")
hist(prevTB$Monocytes,main = "Monocytes")
hist(prevTB$Macrophages.M0,main = "Macrophages.M0")
hist(prevTB$Macrophages.M1,main = "Macrophages.M1")
hist(prevTB$Macrophages.M2,main = "Macrophages.M2")
hist(prevTB$Dendritic.cells.resting,main = "Dendritic.cells.resting")
hist(prevTB$Dendritic.cells.activated,main = "Dendritic.cells.activated")
hist(prevTB$Mast.cells.resting,main = "Mast.cells.resting")
hist(prevTB$Mast.cells.activated,main = "Mast.cells.activated")
hist(prevTB$Eosinophils,main = "Eosinophils")
hist(prevTB$Neutrophils,main = "Neutrophils")

#Reccurent TB immune profile distribution normality test
hist(rectb$B.cells.naive,main = "B.cells.naive")
hist(rectb$B.cells.memory,main = "B.cells.memory")
hist(rectb$Plasma.cells,main = "Plasma.cells")
hist(rectb$T.cells.CD8,main = "T.cells.CD8")
hist(rectb$T.cells.CD4.naive,main = "T.cells.CD4.naive")
hist(rectb$T.cells.CD4.memory.resting,main = "T.cells.CD4.memory.resting")
hist(rectb$T.cells.CD4.memory.activated,main = "T.cells.CD4.memory.activated")
hist(rectb$T.cells.follicular.helper,main = "T.cells.follicular.helper")
hist(rectb$T.cells.regulatory..Tregs.,main = "T.cells.regulatory..Tregs.")
hist(rectb$T.cells.gamma.delta,main = "T.cells.gamma.delta")
hist(rectb$NK.cells.resting,main = "NK.cells.resting")
hist(rectb$NK.cells.activated,main = "NK.cells.activated")
hist(rectb$Monocytes,main = "Monocytes")
hist(rectb$Macrophages.M0,main = "Macrophages.M0")
hist(rectb$Macrophages.M1,main = "Macrophages.M1")
hist(rectb$Macrophages.M2,main = "Macrophages.M2")
hist(rectb$Dendritic.cells.resting,main = "Dendritic.cells.resting")
hist(rectb$Dendritic.cells.activated,main = "Dendritic.cells.activated")
hist(rectb$Mast.cells.resting,main = "Mast.cells.resting")
hist(rectb$Mast.cells.activated,main = "Mast.cells.activated")
hist(rectb$Eosinophils,main = "Eosinophils")
hist(rectb$Neutrophils,main = "Neutrophils")

#LTBI

hist(LTBI$B.cells.naive,main = "B.cells.naive")
hist(LTBI$B.cells.memory,main = "B.cells.memory")
hist(LTBI$Plasma.cells,main = "Plasma.cells")
hist(LTBI$T.cells.CD8,main = "T.cells.CD8")
hist(LTBI$T.cells.CD4.naive,main = "T.cells.CD4.naive")
hist(LTBI$T.cells.CD4.memory.resting,main = "T.cells.CD4.memory.resting")
hist(LTBI$T.cells.CD4.memory.activated,main = "T.cells.CD4.memory.activated")
hist(LTBI$T.cells.follicular.helper,main = "T.cells.follicular.helper")
hist(LTBI$T.cells.regulatory..Tregs.,main = "T.cells.regulatory..Tregs.")
hist(LTBI$T.cells.gamma.delta,main = "T.cells.gamma.delta")
hist(LTBI$NK.cells.resting,main = "NK.cells.resting")
hist(LTBI$NK.cells.activated,main = "NK.cells.activated")
hist(LTBI$Monocytes,main = "Monocytes")
hist(LTBI$Macrophages.M0,main = "Macrophages.M0")
hist(LTBI$Macrophages.M1,main = "Macrophages.M1")
hist(LTBI$Macrophages.M2,main = "Macrophages.M2")
hist(LTBI$Dendritic.cells.resting,main = "Dendritic.cells.resting")
hist(LTBI$Dendritic.cells.activated,main = "Dendritic.cells.activated")
hist(LTBI$Mast.cells.resting,main = "Mast.cells.resting")
hist(LTBI$Mast.cells.activated,main = "Mast.cells.activated")
hist(LTBI$Eosinophils,main = "Eosinophils")
hist(LTBI$Neutrophils,main = "Neutrophils")


#statistical analysis of the immune cell profiles from different patient groups

#Data to be used
attach(prevTB)
attach(rectb)
attach(LTBI)


#computation of Kruskal-Wallis test

colnames(LTBI)

#Multiple pairwise-comparison between gs
pairwise.wilcox.test(LTBI$Mast.cells.activated, LTBI$challeng_group,
                     p.adjust.method = "BH")

#Post-hoc statistical test
aov.ex1 = aov(LTBI$B.cells.naive~LTBI$challeng_group)
summary(aov.ex1,intercept = T)

TukeyHSD(aov.ex1, conf.level=.95)

plot(TukeyHSD(aov(LTBI$Macrophages.M0~LTBI$challeng_group), conf.level=.95))















































