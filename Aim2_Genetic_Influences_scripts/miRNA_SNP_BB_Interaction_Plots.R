# miRNA SNP interaction Plots
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
# Load data ####
load("CandidateGenes_4_2020.RData")

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")
topsnps<-read.csv("Significant_SNPs_GCTA_COJO.csv", stringsAsFactors = F)
topsnpsInfo <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data/TopSNPs.csv")
cluster <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Clustering_miRNA/Figures/miRNA_MEs_nofilter.csv")
names(cluster)[1]<- "shareid"
setwd("/home/clary@mmcf.mehealth.org/Framingham/eQTL Analysis/data")
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

# Format Data ####
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/figure")
data<-pheno_snp_pos %>% unnest(., cols = c(data))
topdata <- data[data$POS %in% topsnps$bp,]
topdata_info <- merge(topdata, topsnpsInfo, by = "POS")
topdata_miRNA <- merge(topdata_info, miRNAdat, by.x = "Indiv", by.y = "shareid")
topdata_cluster <- merge(topdata_miRNA, cluster, by.x = "Indiv", by.y = "shareid")
topdata_all <- topdata_cluster %>% group_by(rsID, POS, Gene.y) %>% nest()
topdata_fem <- topdata_cluster %>% filter(., SEX == "Female") %>% group_by(rsID, POS, Gene.y) %>% nest()
topdata_male <- topdata_cluster %>%  filter(., SEX == "Male") %>% group_by(rsID, POS, Gene.y) %>% nest()

# Plot Function ####
plot_miR_SNP <- function(rsID, miRNA){
  all <- ggplot(data=topdata_all$data[[which(topdata_all$rsID == rsID)]], aes_string(x= "gt_DS", y = miRNA, col = "BB")) + geom_smooth(method = "lm") + 
    ggtitle(paste(topdata_all$rsID[which(topdata_all$rsID == rsID)],topdata_all$Gene.y[which(topdata_all$rsID == rsID)], "Sex-Combined" ,sep = ", " )) + theme(legend.position = "top") +
    scale_color_brewer(palette="Set1")
  female <- ggplot(data=topdata_fem$data[[which(topdata_fem$rsID == rsID)]], aes_string(x= "gt_DS", y = miRNA, col = "BB")) + geom_smooth(method = "lm")  + 
    ggtitle(paste(topdata_fem$rsID[which(topdata_fem$rsID == rsID)],topdata_fem$Gene.y[which(topdata_fem$rsID == rsID)], "Female Only" ,sep = ", " ))+ theme(legend.position = "top") +
    scale_color_brewer(palette="Set1")
  male <- ggplot(data=topdata_male$data[[which(topdata_male$rsID == rsID)]], aes_string(x= "gt_DS", y = miRNA, col = "BB")) + geom_smooth(method = "lm") + 
    ggtitle(paste(topdata_male$rsID[which(topdata_male$rsID == rsID)],topdata_male$Gene.y[which(topdata_male$rsID == rsID)], "Male Only" ,sep = ", " ))+ theme(legend.position = "top") +
    scale_color_brewer(palette="Set1")
  plot <- grid.arrange(all, male, female, nrow = 1)
  ggsave(filename = paste(topdata_all$rsID[which(topdata_all$rsID == rsID)],topdata_all$Gene.y[which(topdata_all$rsID == rsID)], miRNA, "plot.jpg", sep = "_"),plot = plot)
}
# Top SNP-miRNA plots ####  
plot_miR_SNP(rsID = "rs12414657", miRNA = "miR_19a_3p") # ADRB1, miR_19a_3p, female plot 
plot_miR_SNP(rsID = "rs34170507", miRNA = "miR_31_5p") # RANK, miR_31_5p, male plot
plot_miR_SNP(rsID = "rs6567268", miRNA = "let_7g_5p") # RANK, let_7g_5p, male plot 
plot_miR_SNP(rsID = "rs6567268", miRNA = "miR_374a_5p") # RANK, miR_374a_5p, male plot 
plot_miR_SNP(rsID = "rs9533166", miRNA = "miR_19b_3p") # RANK, miR_19b_3p, male plot 
plot_miR_SNP(rsID = "rs9533166", miRNA = "miR_133a") # RANK, miR_133a, male plot
plot_miR_SNP(rsID = "rs9533166", miRNA = "miR_186_5p_a2") # RANK, miR_186_5p_a2, male plot
plot_miR_SNP(rsID = "rs9533166", miRNA = "MEblue") # RANK, MEblue, male plot
plot_miR_SNP(rsID = "rs11124190", miRNA = "miR_17_5p") # HDAC4, miR-17-5p, female plot 
plot_miR_SNP(rsID = "rs11124190", miRNA = "miR_93_5p") # HDAC4, miR-93-5p, female plot 

max(na.omit(cluster[-1]))
