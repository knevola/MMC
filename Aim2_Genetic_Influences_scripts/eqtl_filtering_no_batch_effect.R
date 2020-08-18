# eQTL Filtering and Analysis
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(tidyverse)
options(stringsAsFactors = F)
# Read Files ####
top_snps <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020/Significant_SNPs_GCTA_COJO.csv")
fem_results1 <-read.csv("Model1_miRNA_LMEKIN_eQTL_topSNPs_fem.csv")
male_results1 <- read.csv("Model1_miRNA_LMEKIN_eQTL_topSNPs_male.csv")
fem_results2 <-read.csv("Model2_miRNA_LMEKIN_eQTL_topSNPs_fem.csv")
male_results2 <- read.csv("Model2_miRNA_LMEKIN_eQTL_topSNPs_male.csv")
fem_results3 <-read.csv("Model3_FisherMethod_eQTL_fem.csv")
male_results3 <- read.csv("Model3_FisherMethod_eQTL_male.csv")
fem_results_clust <-read.csv("Cluster_miRNA_LMEKIN_topSNPs_fem.csv")
male_results_clust <- read.csv("Cluster_miRNA_LMEKIN_topSNPs_male.csv")

fem_results1$Model <- "Model_1"
fem_results1$Sex <- "Female"
male_results1$Model <- "Model_1"
male_results1$Sex <- "Male"
fem_results2$Model <- "Model_2"
fem_results2$Sex <- "Female"
male_results2$Model <- "Model_2"
male_results2$Sex <- "Male"
fem_results3$Model <- "Model_3"
fem_results3$Sex <- "Female"
male_results3$Model <- "Model_3"
male_results3$Sex <- "Male"
fem_results_clust$Model <- "Model__clust"
fem_results_clust$Sex <- "Female"
male_results_clust$Model <- "Model__clust"
male_results_clust$Sex <- "Male"

# Format Model_3 results ####
fem_results3 <- fem_results3 %>%  select(., n.x, Value.x, Std.Error.x, z.x, p.value, POS.x, Gene.x, miRNA.x, var.x, Model, Sex)
names(fem_results3) <- names(fem_results2)

male_results3 <- male_results3 %>%  select(., n.x, Value.x, Std.Error.x, z.x, p.value, POS.x, Gene.x, miRNA.x, var.x, Model, Sex)
names(male_results3) <- names(male_results2)

# Combine Results ####
fem_results <- rbind(fem_results1,fem_results2, fem_results3, fem_results_clust)
male_results <- rbind(male_results1, male_results2, male_results3, male_results_clust)

fem_results_clust_sig <- fem_results_clust %>% filter(., p < 0.05) %>% filter(., miRNA %in% c("MEblue", "MEbrown")) %>% filter(., var %in% c("gt_DS:BBYes"))
male_results_clust_sig <- male_results_clust %>% filter(., p < 0.05) %>% filter(., miRNA %in% c("MEblue", "MEbrown")) %>% filter(., var %in% c("gt_DS:BBYes"))

fem_POS <-  c(115803375,239972561,240223080,106736732,60025809)
male_POS <- c(240050108,240112014,84682179,60001153,60026732,43177169)
corr_snps <- read.csv("CorrelationbtwSNPs_4_2020.csv")
osteo_miRs <- read.csv("osteomiRs.csv", header = T)
osteomiRs <- osteo_miRs$X  
osteomiRs <- gsub(pattern = "-",replacement = "_", osteomiRs)

corr_top_snps <- corr_snps %>% filter(., POS1 %in% c(fem_POS, male_POS)| POS2 %in% c(fem_POS, male_POS))
fem_results_SNPs <- fem_results %>%  filter(., var %in% c("gt_DS:BBYes")) %>%filter(., POS %in% fem_POS)
male_results_SNPs <- male_results %>%  filter(., var %in% c("gt_DS:BBYes")) %>% filter(.,POS %in% male_POS)

fem_results_SNPs$FDR <- p.adjust(fem_results_SNPs$p)
male_results_SNPs$FDR <- p.adjust(male_results_SNPs$p)

fem_results_sig <- fem_results_SNPs %>%  filter(., p < 0.05)
male_results_sig <- male_results_SNPs %>% filter(., p < 0.05)

table(male_results_sig$POS,male_results_sig$Gene)
table(fem_results_sig$POS, fem_results_sig$Gene )
