# Fisher P-values Calculation ####
rm(list = ls())
options(stringsAsFactors = F)
library(metaseqR)
library(dplyr)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
fem_results_lin <- read.csv("Model3_linear_miRNA_LMEKIN_eQTL_topSNPs_fem_Batch.csv")
male_results_lin <- read.csv("Model3_linear_miRNA_LMEKIN_eQTL_topSNPs_male_Batch.csv")
fem_results_log <- read.csv("Model3_logistic_miRNA_LMEKIN_eQTL_topSNPs_fem_Batch.csv")
male_results_log <- read.csv("Model3_logistic_miRNA_LMEKIN_eQTL_topSNPs_male_Batch.csv")

fisherpvalue <- function(data1, data2, vars){
  data1 <- data1 %>% filter(.,var == vars)
  data2 <- data2 %>%  filter(.,var == vars)
  data2$miRNA <- gsub("_nas", "", data2$miRNA)
  data1$miRNAPos <- paste(data1$miRNA, data1$POS, sep = ":")
  data2$miRNAPos <- paste(data2$miRNA, data2$POS, sep = ":")
  model1a2<- merge(data1, data2, by= 'miRNAPos')
  row.names(model1a2)<- model1a2$miRNAPos
  keep1 <- c("p.x","p.y") 
  model1a26<- model1a2[names(model1a2) %in% keep1]
  results<-fisher.method(model1a26, p.corr = "BH")
  model3<-merge(model1a2,results, by.x = 0, by.y = 0)
  return(model3)
}

fem_results<-fisherpvalue(data1 = fem_results_lin, data2 = fem_results_log, vars = "gt_DS:BBYes")
male_results<-fisherpvalue(data1 = male_results_lin, data2 = male_results_log, "gt_DS:BBYes")

min(fem_results$p.value)
min(fem_results$p.adj)
min(male_results$p.value)
min(male_results$p.adj)

fem_results_sig <-fem_results %>% filter(., p.value < 0.05)
male_results_sig <- male_results %>% filter(., p.value < 0.05)

write.csv(fem_results, "Model3_FisherMethod_eQTL_fem_Batch.csv", quote = F, row.names = F)
write.csv(male_results, "Model3_FisherMethod_eQTL_male_Batch.csv", quote = F, row.names = F)
