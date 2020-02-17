library(scales)
library(dplyr)
library(caret)
library(devtools)
library(ggplot2)
library(tidyverse)
library(tseries)
library(reshape2)
library(imputeTS)
library(psych)
library(Rtsne)
library(GPArotation)
library(cluster)
library(nonlinearTseries)

dataset <- read.csv("Data/ESMdata.csv", header = TRUE)
dataset$dayno <- dataset$dayno[order(dataset$dayno)]
dataset$dayno <- as.numeric(factor(dataset$dayno)) # re factor 'dayno' due to gaps in 'dayno'

newDS <- na_interpolation(dataset)# Impute the missing values with na.interpolation
df <- data.frame(newDS)
options(stringsAsFactors = FALSE)
df <- df %>% group_by(dayno)
#################################

#detecting, removing influential outliers
outlierdf <- (df[c(11:30)])
#outlier viz cook dist
cd <- cooks.distance(lm(outlierdf))
cddf <- as.data.frame(cd)
cddf$index <- seq.int(nrow(cddf))

ggplot(cddf, aes(x = index, y = cd)) +
   geom_point(color="red") +
   geom_hline(yintercept = 4*mean(cd), color="blue") 

#removed influentials
influential <- as.numeric(names(cd)[(cd > 4*mean(cd, na.rm=T))]) 
infRmDf <- df[-influential, ]

#reduce n dimenssions to 2 for ggplot (there are better packages than ggplot for PCA/EFA etc)
tsne.pca = Rtsne(as.matrix(infRmDf[c(11:30)]), scale = TRUE, check_duplicates=FALSE, pca=TRUE, perplexity=100, theta=0.5, dims=2)

tsne_V1V2 = as.data.frame(tsne.pca$Y)

tsne_original <- tsne_V1V2

#kmeans - optimal groups -> 3 (Eigenvalue >1)
fit_cluster_kmeans <-  kmeans(scale(tsne_V1V2), 4)
tsne_original$cl_kmeans <-  factor(fit_cluster_kmeans$cluster)

tsne_V1V2$cluster <- tsne_original$cl_kmeans

ggplot(data=tsne_V1V2, aes(x=V1, y=V2, color=tsne_original$cl_kmeans)) + 
  geom_point() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mental State Clusters\n", x = "Dim 1", y = "Dim 2", color = "Mental State Clusters\n")  

fit_cluster_kmeans$size

####################
#Corrolation matrix 20 mood types
fa_cor <- cor(infRmDf[c(11:30)], method = "spearman")
fa_cor <- sqrt(fa_cor^2) #try for normality
fa_cor = as.data.frame(t(fa_cor))
fa_cor$groups <- NA
fa_cor$RowNumber <- seq.int(nrow(fa_cor))

####################
#Assumption Test for FA/PCA  
#cortest.bartlett(infRmDf[c(11:30)] , n=1476) #(Snedecor et al 1989)
#out -> Bartlett's test is significant, x2(190) = 23514.21, p < .001 - (Cerny & Kaiser 1977, pp.43-47)
# and therefore factor analysis is appropriate.
#KMO(r = infRmDf[, c(11:30)]) #Overall MSA =  0.94 
#fa.parallel(infRmDf[,c(11:30)],fm="pa") #quasi 4
#alpha(df[c(11:30)], check.keys=TRUE) #Cronbach alpha, >0.94 - (Yurdugül, H, 2008)

# negitive loadings - high score on the factor is associated with a lower score on the item.

#residuals view 
#residuals <- factor.residuals(fa_cor[c(1:20)], faGroups$loadings)
#residuals <- as.matrix(residuals[upper.tri(residuals)])
#large.resid <- abs(residuals) > 0.05
#sqrt(mean(residuals^2))
#ggplot(residuals, aes(x=residuals))+ 
#geom_histogram(color="black", fill="white", binwidth=.1)

Dist <- daisy(fa_cor[c(1:20)], metric = c("gower"), type=list())
Dist <- as.matrix(Dist)
heatmap(Dist, Rowv=TRUE, Colv="Rowv", symm = TRUE)

#PCA (rotated) to understand component groups
#################################
#Eigenvalue >1 (No.3, similar to parallel test)
eign.df <- psych::principal(infRmDf[c(11:30)], nfactors=20)
ggplot(, aes(x=fa_cor$RowNumber, y=eign.df$values)) + 
  geom_line() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Eigenvalues > 1\n", x = "Factors", y = "Eigenvalues") +
  geom_hline(yintercept = 1, color="blue")  +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  geom_point()

fa.pca <- psych::principal(infRmDf[c(11:30)], nfactors=4, rotate="oblimin", scores = TRUE)
fa.groups <- print.psych(fa.pca , cut = 0.43, sort = TRUE)
#fa.diagram(fa.pca)

#Assign groups based on PCA routed / further PCA analaysys not done here
#fa_cor[fa_cor$RowNumber %in% c(16, 20, 12, 9, 17, 11, 7, 4, 19), "groups"] <- "Mood Positive"
#fa_cor[fa_cor$RowNumber %in% c(3, 14, 13, 1), "groups"] <- "Mood Uncertain"
#fa_cor[fa_cor$RowNumber %in% c(15, 8, 18, 2), "groups"] <- "Mood Negitive"
#fa_cor[fa_cor$RowNumber %in% c(6, 10, 5), "groups"] <- "Mood Worried"

#mood status summarise
#########################
moods <- infRmDf
# negative affect
Mood_Positive <- c(26, 30,22, 19, 27, 21, 17, 14, 29)
moods$mp <-  moods %>% select(Mood_Positive) %>% rowMeans(.,na.rm = TRUE)

Mood_Uncertain <- c(13, 24, 23, 11)
moods$mu <- moods %>% select(Mood_Uncertain) %>% rowMeans(.,na.rm = TRUE)

Mood_Negitive <- c(25, 28, 28, 12)
moods$mn <- moods %>% select(Mood_Negitive) %>% rowMeans(.,na.rm = TRUE)

Mood_Worried <- c(16, 20, 15)
moods$mw <- moods %>% select(Mood_Worried) %>% rowMeans(.,na.rm = TRUE)

moods$moodsums <- rowSums(moods[c("mp", "mu", "mn", "mw")])

moods$groups <- as.numeric(tsne_original$cl_kmeans)

#Group day numbers and mood means
moods.dfa <- moods %>% group_by(dayno) %>% summarise(moodsums = mean(moodsums), groups = floor(mean(groups)))
moods.dfa$moodsums_dfa = NA

# determine DFA, in a moving window of 20 days
window <- 21
for (i in seq(window, max(moods$dayno), 1)) {

  # get the sliding window data
  sliding.window <- subset(moods, dayno > (i - window) & dayno <= i)
    
  # moodsums
  dfa.anal <- dfa(time.series = sliding.window$moodsums,
                      npoints = 30,
                      window.size.range = c(10, nrow(sliding.window)),
                      do.plot = FALSE)
  
  fgn.est <- estimate(dfa.anal) 

  moods.dfa$moodsums_dfa[moods.dfa$dayno == i] <- fgn.est
  }
fit_cluster_kmeans$size
#assessment periods
rects <- data.frame(xstart = c(1, 29,43,99,156), xend = c(28,42,98,155,238)) #(Kossakowski et al., 2017)

ggplot(data = na.omit(moods.dfa$moods.dfa)) + 
      geom_point() +
      geom_rect(data=rects, aes( xmin=xstart, xmax=xend, ymin = -Inf, ymax = Inf), alpha =0.1) + 
      geom_line(aes(x=moods.dfa$dayno, y=moods.dfa$moodsums_dfa, colour= factor(moods.dfa$groups),  group=1 ), ) +
      ylab("DFA (21-day window)") + 
      xlim(c(0, 250)) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Mental States over Test Period with phases (shadowed) \n", x = "Test Day Number", y = "21 day Window" )  +
      scale_color_discrete(name = "Mental State Groups",  labels=c("Mood Positive", "Mood Uncertain", "Mood Worried", "Mood Negitive"))    +  
      annotate("text", x = 15, y = 2.08, label = "Baseline") +  
      annotate("text", x = 36, y = 2.06, label = "No Change") +
      annotate("text", x = 70, y = 2.08, label = "Reduce AD") +
      annotate("text", x = 125, y = 2.08, label = "Post test") +
      annotate("text", x = 195, y = 2.08, label = "Follow up")
  
  
   