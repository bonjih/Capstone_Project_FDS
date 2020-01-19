library(recipes)
library(dplyr)
library(caret)
library(devtools)
library(ggplot2)
library(ggbiplot)
library(tidyverse)
library(tseries)
library(naniar)
library(reshape2)
library(ggfortify)
library(imputeTS)
library(lubridate)
library(magrittr)
library(data.table)

#rm(list = ls()) 

dataset <- read.csv("data/ESMdata.csv", header = TRUE)
newDS <- na.interpolation(dataset)# Impute the missing values with  na.interpolation
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

#PCA
mood.pca <- df[c(11:30)]
dfCors <- cor(mood.pca[ ,1:20], method="spearman") # convert to a value between 0 & 1
df.pca <- prcomp(dfCors [,c(1:20)], center = TRUE, scale. = TRUE) 
summary(df.pca)
ggbiplot(df.pca, labels = rownames(dfCors))

