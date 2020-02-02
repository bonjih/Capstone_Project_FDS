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


rm(list = ls()) 

dataset <- read.csv("Data/ESMdata.csv", header = TRUE)

newDS <- na_interpolation(dataset)# Impute the missing values with na.interpolation
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

#################################

#detecting, removing outliers
outlierdf <- (df[c(11:30)])
#outlier viz cook dist
cd <- cooks.distance(lm(outlierdf))
plot(cd,pch=1, cex = 1)
abline(h = 4*mean(cd,na.rm = T), col = "blue")
text(x=1:length(cd)+2, y=cd, labels=ifelse(cd>4*mean(cd, na.rm=T),names(cd),""), col="blue")
#remove
influential <- as.numeric(names(cd)[(cd > 4*mean(cd, na.rm=T))]) 
infRmDf <- df[-c(1, 21, 41), ]

#EFA to fine factors

efa_cor <- cor(infRmDf[c(11:30)], method = 'spearman')
efa_cor = as.data.frame(t(efa_cor))
efa_cor$groups <- NA
efa_cor$RowNumber <- seq.int(nrow(efa_cor))

####################
#Assumption Test  

#cortest.bartlett(df[c(11:30)] , n=1476) (Snedecor et al 1989)
#out -> Bartlett's test is highly significant, x2(190) = 23547, p < .001 - (Cerny & Kaiser 1977, pp.43-47)
# and therefore factor analysis is appropriate.
#KMO(r = df[, c(11:30)]) #Overall MSA =  0.95 
#fa.parallel(df[,c(11:30)],fm="pa") #quasi 3
#alpha(df[c(11:30)], check.keys=TRUE) 
#Cronbach alpha, >0.94 - (Yurdugül, H, 2008)

####################

efaGroups <- principal(efa_cor[c(1:20)], nfactors=3, rotate="oblimin", scores = TRUE)
efa.groups <- print.psych(efaGroups , cut = 0.4, sort = TRUE, scores = TRUE)

#Assign groups based on efaGroups
efa_cor[efa_cor$RowNumber %in% c(16, 12, 17, 20, 9, 4, 11, 7, 19, 18), "groups"] <- "Mood Positive"
efa_cor[efa_cor$RowNumber %in% c(5, 6, 10, 2, 15, 8), "groups"] <- "Mood Negitive"
efa_cor[efa_cor$RowNumber %in% c(3, 14, 13,1), "groups"] <- "Mood Uncertain"

#use EFA and kmeans confirm number of gruops
#kmeans - optimal groups -3
DataStandardised <- scale(infRmDf[11:30])
DataStandardised <-data.matrix(DataStandardised)
groups <- kmeans(DataStandardised, 3)
clusplot(DataStandardised, groups$cluster)
str(groups)

fviz_cluster(groups, data = DataStandardised )  +
  ggtitle("Kmeans Mental State CLusters") +
  theme(plot.title = element_text(hjust = 0.5))   
 
#################################

#PCA setup
#mood <- df[c(2,3, 5,11:30)]
#mood <- mood %>% group_by(date) %>% summarise_at(vars(phase:se_handle), mean) 
#moodDaily <- mood[c(4:22)] 

#########################

moodefa_cor.pca <- prcomp(efa_cor[c(1:20)], center = TRUE) # includes data standardisation

#ggbiplot::ggbiplot(moodefa_cor.pca, ellipse=TRUE, labels=rownames(efa_cor), groups=efa_cor$groups)  

#moodefa_cor.pca$rotation

#PCA Eigenvalue
screeplot(moodefa_cor.pca, type = "l", npcs = 5, main = "Screeplot of the first 5 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

#par(mfrow=c(2,2))
 
#summary(moodefa_cor.pca)

#PCA Plot
fviz_pca_biplot(moodefa_cor.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = efa_cor$groups,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Mental State") +
  ggtitle("PCA plot - Participant's Mental State Over Test Period") +
  theme(plot.title = element_text(hjust = 0.5))








