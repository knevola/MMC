# Bar graph treatment for hypertension
rm(list = ls())

library(tidyverse)
library(gt)
library(dplyr)
library(gtsummary)
library(ggplot2)
library(emmeans)
library(reshape2)
library(matrixStats)
library(haven)

# Read in data
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')

miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

miRNA_pheno <- merge(pheno, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$miR19a <- miRNA_pheno$miR_19a_3p
miRNA_pheno$miR186 <- miRNA_pheno$miR_186_5p_a2
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "Yes" & miRNA_pheno$BB == "Yes"] <-"BB Hypertension Treatment"
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "Yes" & miRNA_pheno$BB == "No"] <-"Non-BB Hypertension Treatment"
miRNA_pheno$HRX_BB[miRNA_pheno$HRX8 == "No"] <-"Non-Hypertension Treatment"
miRNA_pheno$HRX_BB <- as.factor(miRNA_pheno$HRX_BB)
miRNA_pheno$HRX_BB <- relevel(x = miRNA_pheno$HRX_BB, ref = "Non-Hypertension Treatment")

