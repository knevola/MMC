#FN BMD ~ SNP Interaction Plots
rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
# Load Files ####
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
load("CandidateGenes_4_2020.RData")

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")
topsnps<-read.csv("Significant_SNPs_GCTA_COJO.csv", stringsAsFactors = F)
topsnpsInfo <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data/TopSNPs.csv")

# Format Data ####
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/figure")
data<-pheno_snp_pos %>% unnest(., cols = c(data))
topdata <- data[data$POS %in% topsnps$bp,]
topdata_info <- merge(topdata, topsnpsInfo, by = "POS")
topdata_all <- topdata_info %>% group_by(rsID, POS, Gene.y) %>% nest()
topdata_fem <- topdata_info %>% filter(., SEX == "Female") %>% group_by(rsID, POS, Gene.y) %>% nest()
topdata_male <- topdata_info %>%  filter(., SEX == "Male") %>% group_by(rsID, POS, Gene.y) %>% nest()

# Figures
for (i in 1:length(topdata_all$rsID)){
  all <- ggplot(data=topdata_all$data[[i]], aes(x= gt_DS, y = f8cbnbmd, col = BB)) + geom_smooth(method = "lm") + ylab("FN BMD") + 
    ggtitle(paste(topdata_all$rsID[i],topdata_all$Gene.y[i], "Sex-Combined" ,sep = ", " )) + theme(legend.position = "top")
  female <- ggplot(data=topdata_fem$data[[i]], aes(x= gt_DS, y = f8cbnbmd, col = BB)) + geom_smooth(method = "lm") + ylab("FN BMD") + 
    ggtitle(paste(topdata_fem$rsID[i],topdata_fem$Gene.y[i], "Female Only" ,sep = ", " ))+ theme(legend.position = "top")
  male <- ggplot(data=topdata_male$data[[i]], aes(x= gt_DS, y = f8cbnbmd, col = BB)) + geom_smooth(method = "lm") + ylab("FN BMD") + 
    ggtitle(paste(topdata_male$rsID[i],topdata_male$Gene.y[i], "Male Only" ,sep = ", " ))+ theme(legend.position = "top")
  plot <- grid.arrange(all, male, female, nrow = 1)
  ggsave(filename = paste(topdata_all$rsID[i],topdata_all$Gene.y[i], "FN_BMD_plot.jpg", sep = "_"),plot = plot)
}

hist(topdata_all$data[[7]]$gt_DS,breaks = 10)
x<-as.data.frame(table(topdata_all$data[[7]]$gt_DS))

ggplot(data = topdata_all$data[[7]],aes(x = gt_DS)) + geom_histogram(bins = 39)

ggplot(data = x,aes(x = Var1, y = Freq)) + geom_bar(stat="identity")
