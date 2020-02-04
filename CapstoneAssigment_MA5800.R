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
library(cluster)


rm(list = ls()) 

dataset <- read.csv("Data/ESMdata.csv", header = TRUE)

newDS <- na_interpolation(dataset)# Impute the missing values with na.interpolation
df <- data.frame(newDS)
options(stringsAsFactors = FALSE)

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

####################
#Corrolation matrix

fa_cor <- cor(infRmDf[c(11:30)], method = "spearman")
fa_cor <- sqrt(fa_cor^2) #try for normality
fa_cor = as.data.frame(t(fa_cor))
fa_cor$groups <- NA
fa_cor$RowNumber <- seq.int(nrow(fa_cor))

####################
#Assumption Test  
#cortest.bartlett(df[c(11:30)] , n=1476) #(Snedecor et al 1989)
#out -> Bartlett's test is significant, x2(190) = 23547, p < .001 - (Cerny & Kaiser 1977, pp.43-47)
# and therefore factor analysis is appropriate.
#KMO(r = df[, c(11:30)]) #Overall MSA =  0.95 
#fa.parallel(df[,c(11:30)],fm="pa") #quasi 3
#alpha(df[c(11:30)], check.keys=TRUE) 
#Cronbach alpha, >0.94 - (Yurdugül, H, 2008)
#det(fa_cor)
#which(fa_cor > 0.8)
####################

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

faGroups <- principal(df[c(11:30)], nfactors=3, rotate="oblimin", scores = TRUE)
fa.groups <- print.psych(faGroups , cut = 0.4, sort = TRUE)
fa.diagram(faGroups)

#residuals view 
residuals <- factor.residuals(fa_cor[c(1:20)], faGroups$loadings)
residuals <- as.matrix(residuals[upper.tri(residuals)])
large.resid <- abs(residuals) > 0.05
sum(large.resid)/nrow(residuals)
sqrt(mean(residuals^2))
hist(residuals)

#Assign groups based on efaGroups
fa_cor[fa_cor$RowNumber %in% c(16, 12, 20, 9, 17, 11, 4, 7, 19), "groups"] <- "Mood Positive"
fa_cor[fa_cor$RowNumber %in% c(5, 10, 2, 15, 6, 8, 18), "groups"] <- "Mood Negitive"
fa_cor[fa_cor$RowNumber %in% c(3, 14, 13,1), "groups"] <- "Mood Uncertain"

#Eigenvalue >1
plot(faGroups$values, type = "b")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)


#########################

#PCA setup
#mood <- df[c(2,3, 5,11:30)]
#mood <- mood %>% group_by(date) %>% summarise_at(vars(phase:se_handle), mean) 
#moodDaily <- mood[c(4:22)] 

moodfa_cor.pca <- prcomp(fa_cor[c(1:20)], center = TRUE) # for fviz_pca_biplot / includes data standardisation
moodfa_cor.pca$x
#ggbiplot::ggbiplot(moodfa_cor.pca, ellipse=TRUE, labels=rownames(fa_cor), groups=fa_cor$groups)  

#moodefa_cor.pca$rotation

#par(mfrow=c(2,2))
 
#summary(moodfa_cor.pca)

#PCA Plot
fviz_pca_biplot(moodfa_cor.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = fa_cor$groups,
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Mental State") +
  ggtitle("PCA plot - Participant's Mental State Over Test Period") +
  theme(plot.title = element_text(hjust = 0.5))


#########################

#split dataset into phases,


Mood.Positive <- mean(c(16, 12, 17, 20, 9, 4, 11, 7, 19, 18))
Mood.Negitive <- mean(c(5, 6, 10, 2, 15, 8))
Mood.Uncertain <- mean(c(3, 14, 13,1))



