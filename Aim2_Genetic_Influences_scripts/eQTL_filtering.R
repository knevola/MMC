# eQTL Filtering and Analysis
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(tidyverse)
options(stringsAsFactors = F)
# Read Files ####
top_snps <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020/Significant_SNPs_GCTA_COJO.csv")
fem_results1 <-read.csv("Model1_miRNA_LMEKIN_eQTL_topSNPs_fem_Batch.csv")
male_results1 <- read.csv("Model1_miRNA_LMEKIN_eQTL_topSNPs_male_Batch.csv")
fem_results2 <-read.csv("Model2_miRNA_LMEKIN_eQTL_topSNPs_fem_Batch.csv")
male_results2 <- read.csv("Model2_miRNA_LMEKIN_eQTL_topSNPs_male_Batch.csv")
fem_results3 <-read.csv("Model3_FisherMethod_eQTL_fem_Batch.csv")
male_results3 <- read.csv("Model3_FisherMethod_eQTL_male_Batch.csv")
fem_results_clust <-read.csv("Cluster_miRNA_LMEKIN_topSNPs_fem_Batch.csv")
male_results_clust <- read.csv("Cluster_miRNA_LMEKIN_topSNPs_male_Batch.csv")

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

# Combining Results to see directional trends
miRNA_fem_sig <- data.frame(miRNA = fem_results_sig$miRNA, POS = fem_results_sig$POS, Gene = fem_results_sig$Gene,Model = fem_results_sig$Model,
                            Estimate = fem_results_sig$Value, SE = fem_results_sig$Std.Error, miRNA_p = fem_results_sig$p, Sex = "females")
miRNA_male_sig <- data.frame(miRNA = male_results_sig$miRNA, POS = male_results_sig$POS, Gene = male_results_sig$Gene, Model = male_results_sig$Model,
                            Estimate = male_results_sig$Value, SE = male_results_sig$Std.Error, miRNA_p = male_results_sig$p, Sex = "males")
miRNA_sig <- rbind(miRNA_fem_sig, miRNA_male_sig)
miRNA_sig$osteomiR <- FALSE
miRNA_sig$osteomiR[miRNA_sig$miRNA %in% osteomiRs]<- TRUE

miRNA_sig_POS <- merge(miRNA_sig, top_snps, by.x = "POS", by.y = "bp")
miRNA_sig_POS$Gene <- gsub("TNFS11","TNFSF11",miRNA_sig_POS$Gene)

# Does the miRNA target the Gene?
library(multiMiR)
miRs<- get_multimir(org = "hsa",target = unique(miRNA_sig_POS$Gene))
miRs_target<- unique(miRs@data %>% dplyr::select(., mature_mirna_id, target_symbol))
miRs_target$mature_mirna_id <- gsub(pattern = "hsa-", replacement = "", miRs_target$mature_mirna_id)
miRs_target$mature_mirna_id <- gsub(pattern = "-", replacement = "_", miRs_target$mature_mirna_id)
miRs_target$Combination <- paste(miRs_target$mature_mirna_id, ":", miRs_target$target_symbol)

miRNA_sig_POS$Combination <- paste(miRNA_sig_POS$miRNA, ":", miRNA_sig_POS$Gene)
miRNA_sig_POS$targetted <- FALSE
miRNA_sig_POS$targetted[miRNA_sig_POS$Combination %in% miRs_target$Combination]<- TRUE
table(miRNA_sig_POS$targetted, miRNA_sig_POS$Gene)

# Filtering SNPs to determine SNPs for validation ####
top_snps$bp[!(top_snps$bp %in% unique(miRNA_sig_POS$POS))]
table(miRNA_sig_POS$POS,miRNA_sig_POS$Sex)
miRNA_sig_POS_01 <- miRNA_sig_POS %>% filter(., pJ < 0.01) 
table(miRNA_sig_POS_01$POS)
miRNA_sig_POS_01_m <- miRNA_sig_POS %>% filter(., miRNA_p < 0.01) 
table(miRNA_sig_POS_01_m$POS, miRNA_sig_POS_01_m$Gene)

# Associated with SNPs of Interest ####
miRNA_sig_POS$miR_interest <- FALSE
miRNA_sig_POS$miR_interest[miRNA_sig_POS$osteomiR == T] <- T
miRNA_sig_POS$miR_interest[miRNA_sig_POS$miRNA %in% c("miR_19a_3p", "miR_186_5p_a2")] <- T
table(miRNA_sig_POS$POS, miRNA_sig_POS$miR_interest)

table(miRNA_sig_POS$POS, miRNA_sig_POS$targetted)
table(miRNA_sig_POS$POS)
table(miRNA_sig_POS_01_m$POS)

# miRNA significant in top SNPs for validation #####
rs12414657_miRNA <- miRNA_sig_POS %>% filter(., POS == 115803375) %>% filter(., targetted == T | miR_interest == T)
rs34170507_miRNA <- miRNA_sig_POS %>% filter(., POS == 60001153) %>% filter(., targetted == T | miR_interest == T)
rs6567268_miRNA <- miRNA_sig_POS %>% filter(., POS == 60026732) %>% filter(., targetted == T | miR_interest == T)
fem_snp_HDAC4 <-miRNA_sig_POS %>%  filter(.,  POS == 240223080) %>% filter(., targetted == T | miR_interest == T)
top_miRs<-rbind(rs12414657_miRNA,rs34170507_miRNA, rs6567268_miRNA, fem_snp_HDAC4)
table(miRNA_sig_POS$POS, miRNA_sig_POS$miR_interest)
table(miRNA_sig_POS_01_m$POS)
write.csv(top_miRs, "Sig_Assoc_MiRNA_top_SNPs.csv", quote = F, row.names = F)
