# Top SNPs eQTL miRNA enrichment analysis
rm(list = ls())
library(tidyverse)
library(multiMiR)
library(limma)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

# Read files
sig_mirs<-read.csv("Female_Top2_SNPs_eQTL_sig.csv")

# Format miRNAs to be read by multimiR
sig_mirs$miRNA <- gsub(pattern = "miR", replacement = "hsa-miR",x = sig_mirs$miRNA)
sig_mirs$miRNA <- gsub(pattern = "_", replacement = "-", x=sig_mirs$miRNA)

# Separate into SNPs
ADRB1_mirs <- sig_mirs %>% filter(., sig_mirs$Gene == "ADRB1")
HDAC4_mirs <- sig_mirs %>% filter(., sig_mirs$Gene == "HDAC4")

# Determine targets
ADRB1_miR_targets<- get_multimir(org = "hsa",mirna = ADRB1_mirs$miRNA)
ADRB1_miR_targets <- ADRB1_miR_targets@data
ADRB1_miR_targets <- unique(ADRB1_miR_targets[3:5])

HDAC4_miR_targets<- get_multimir(org = "hsa",mirna = HDAC4_mirs$miRNA)
HDAC4_miR_targets <- HDAC4_miR_targets@data
HDAC4_miR_targets <- unique(HDAC4_miR_targets[3:5])

# Filter for targetted by more that 1 miRNA
ADRB1_target_freq<-as.data.frame(table(ADRB1_miR_targets$target_entrez))
ADRB1_target_freq_3 <- ADRB1_target_freq[(ADRB1_target_freq$Freq > 2) & (ADRB1_target_freq$Freq < 10),]
ADRB1_target_freq_3$Var1 <- as.character(ADRB1_target_freq_2$Var1)

HDAC4_target_freq<-as.data.frame(table(HDAC4_miR_targets$target_entrez))
HDAC4_target_freq_3 <- HDAC4_target_freq[(HDAC4_target_freq$Freq > 2) & (HDAC4_target_freq$Freq < 10),]
HDAC4_target_freq_3$Var1 <- as.character(HDAC4_target_freq_2$Var1)

# Enrichment analysis
ADRB1_GO <- goana(de = unique(ADRB1_target_freq_3$Var1))
ADRB1_GO_sig <- ADRB1_GO %>% filter(., P.DE < 0.05)
ADRB1_KEGG <- kegga(de = ADRB1_target_freq_3$Var1)
ADRB1_KEGG_sig <- ADRB1_KEGG %>% filter(., P.DE < 0.05)

HDAC4_GO <- goana(de = unique(HDAC4_target_freq_3$Var1))
HDAC4_GO_sig <- HDAC4_GO %>% filter(., P.DE < 0.05)
HDAC4_KEGG <- kegga(de = HDAC4_target_freq_3$Var1)
HDAC4_KEGG_sig <- HDAC4_KEGG %>% filter(., P.DE < 0.05)
