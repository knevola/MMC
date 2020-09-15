rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(tidyverse)

load("CandidateGenes_4_2020.RData")
data<-pheno_snp_pos %>% unnest(., cols = c(data))
correlation <- NULL

`%notin%` <- Negate(`%in%`)

data_nested<-data %>% group_by(., Gene, POS) %>% nest() 
for (g in  unique(data_nested$Gene)){
  print(g)
  data1 <- data_nested %>% filter(., Gene == g)
  for (i in 1:length(data1$POS)){
    print(i)
    for(j in i:length(data1$POS)){
      if (data1$POS[i] != data1$POS[j]){
        c <- cor(data1$data[[i]]$gt_DS,data1$data[[j]]$gt_DS, method = "spearman")
        if (c > 0.9){
          corrs <- data.frame(Gene = g, POS1 = data1$POS[i], POS2 = data1$POS[j], corr = c)
          correlation <- rbind(correlation, corrs)
        }
      }
    }
  }
}

write.csv(correlation, "CorrelationbtwSNPs_4_2020.csv", quote = F, row.names = F)

correlation <- read.csv("CorrelationbtwSNPs_4_2020.csv")

corr1 <- correlation[correlation$corr == 1,]
length(unique(corr1$POS1)) #8
length(unique(corr1$POS2)) #8

corr99 <- correlation[correlation$corr > .99,]
length(unique(corr99$POS1)) #421
length(unique(corr99$POS2)) #422

corr98 <- correlation[correlation$corr > .98,]
length(unique(corr98$POS1)) #632
length(unique(corr98$POS2)) #631

corr97 <- correlation[correlation$corr > .97,]
length(unique(corr97$POS1)) #748
length(unique(corr97$POS2)) #748

corr96 <- correlation[correlation$corr > .96,]
length(unique(corr96$POS1)) #818
length(unique(corr96$POS2)) #821

corr95 <- correlation[correlation$corr > .95,]
length(unique(corr95$POS1)) #874
length(unique(corr95$POS2)) #880


data96 <- data %>% filter(., POS %notin% corr96$POS2) %>% group_by(., POS, Gene) %>% nest()
save(data96, file = "allSNPs_cleaned_corr96.RData")
table(data96$Gene)


