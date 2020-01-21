library(ggbiplot)
library(scales)
library(dplyr)
library(recipes)
library(caret)
library(devtools)
library(ggplot2)
library(tidyverse)
library(tseries)
library(naniar)
library(reshape2)
library(ggfortify)
library(imputeTS)
library(lubridate)
library(magrittr)
library(rlang)
library(Hmisc)
library(latticeExtra)
library(factoextra)
library(lattice)
library(FactoMineR)
library(psych)
library(GPArotation)

#rm(list = ls()) 

dataset <- read.csv("data/ESMdata.csv", header = TRUE)
newDS <- na_interpolation(dataset)# Impute the missing values with  na.interpolation
df <- data.frame(newDS)
options(stringsAsFactors = FALSE)

#Time series for 'dep' results
dfts <- df %>%
  mutate(date = as.Date(date, format="%d/%m/%Y")) %>%
  complete(date = full_seq(date, period = 1)) 

ggplot(dfts, aes(x = date, y = dep)) +
  geom_line() +
  scale_x_date(name = "Date (239 days)", 
               date_labels = "%d-%m-%Y", 
               date_breaks = "2.5 weeks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#PCA setup
mood <- df[c(2,3, 5,11:30)]
mood <- mood %>% group_by(date) %>% summarise_at(vars(phase:se_handle), mean) 
moodDaily <- mood[c(4:22)]

rownames(mood) = make.names(mood$phase, unique=TRUE )

#moodCors <- cor(moodDaily, method="spearman") # convert to a value between 0 & 1
moodDaily.pca <- prcomp(moodDaily, center = TRUE, scale. = TRUE)

moodDaily.pca$rotation

#PCA Eigenvalue
screeplot(moodDaily.pca, type = "l", npcs = 5, main = "Screeplot of the first 5 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

plot(moodDaily.pca$x[,1],moodDaily.pca$x[,2], xlab="PC1 (60.6%)", ylab = "PC2 (11.4%)", main = "PC1 / PC2 - plot")
 
#summary(moodDaily.pca)

#PCA Plot
#ggbiplot::ggbiplot(moodDaily.pca, labels = rownames(mood), 
                   #circle = TRUE, var.scale = 1 ,
                   #obs.scale = 1, varname.size = 4) 
                   #ggtitle("PCA of Patient moods - 238 obs")  
              
biplot(moodDaily.pca, col=c("black","red"), cex=c(0.7,0.8))

fviz_pca_biplot(moodDaily.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = as.factor(mood$phase),
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Test Phase") +
  ggtitle("PCA plot Patient mood v test Phase") +
  theme(plot.title = element_text(hjust = 0.5))








