# Table 1: Characteristic of Study sample #
rm(list = ls())
# Install packages if you have not previously installed:
# install.packages("tidyverse")
#devtools::install_github("rstudio/gt")
#install.packages("remotes")
#remotes::install_github("ddsjoberg/gtsummary")

library(tidyverse)
library(gt)
library(dplyr)
library(gtsummary)
library(ggplot2)
library(emmeans)
library(haven)

# Read in data
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('GEN3Pheno.csv')
miRNA_pheno <- merge(pheno, miRNAdat, by.x = "shareid", by.y = "shareid")
pheno1 <- miRNA_pheno[2:31]
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

# Merge pheno with miRNA ####
quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)
Colsums <- data.frame(variables = names(miRNA_pheno), NAs = colSums(is.na(miRNA_pheno))/nrow(miRNA_pheno))
Colsums$NAs[1:30]<- NA 
drop <- c('cvdpair', 'casecontrol', 'idtype.x', "idtype.y")
miRNA_pheno <- miRNA_pheno[, !(names(miRNA_pheno) %in% drop)]

fnmodel <- lm(f2nbmd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
ftomodel<-lm(f2tobmd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
ftrmodel <- lm(f2trbmd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
s2model <- lm(s2l2bd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
s3model <- lm(s2l3bd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
s4model <- lm(s2l4bd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)
s24model <- lm(s2l24bd~BB + AGE2 + SEX + HGT2 + WGT2, data = pheno1)

summary(fnmodel)
summary(ftomodel)
summary(ftrmodel)
summary(s2model)
summary(s3model)
summary(s4model)
summary(s24model)

#miR_19a_3p and miR_186_5p_a2 validation
mlr19a_3bmd <- lm(miR_19a_3p~s2l24bd + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
mlr19a_3bmd_f <- lm(miR_19a_3p~f2tobmd + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
mlr19a_3bb <- lm(miR_19a_3p~BB + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
summary(mlr19a_3bmd)
summary(mlr19a_3bmd_f)
summary(mlr19a_3bb)

mlr186_3bmd <- lm(miR_186_5p_a2~s2l24bd + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
mlr186_3bmd_f <- lm(miR_186_5p_a2~f2tobmd + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
mlr186_3bb <- lm(miR_186_5p_a2~BB + AGE2 + SEX + HGT2 + WGT2 + rankqual + rankcon + rank260, data = miRNA_pheno)
summary(mlr186_3bmd)
summary(mlr186_3bmd_f)
summary(mlr186_3bb)
