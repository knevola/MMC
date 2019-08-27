#7/8 #
# Setup ####
rm(list = ls())
library(tidyverse)
library(dplyr)
library(VennDiagram)
library(metaseqR)
library(ggplot2)
library(pheatmap)
library(haven)
library(lmerTest)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

# Merge pheno with miRNA ####
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)

miRNA_pheno$miR_19a_3p
miRNA_pheno$miR_186_5p_a2

#Distribution of Technical Variables
hist(miRNA_pheno$concentration)
miRNA_pheno$log_con <- log(miRNA_pheno$concentration)
hist(miRNA_pheno$log_con)# More normal still skewed

hist(miRNA_pheno$`_260_280`)
miRNA_pheno$new260 <- ifelse((miRNA_pheno$`_260_280` < 1 | miRNA_pheno$`_260_280` >3), NA, miRNA_pheno$`_260_280`)
hist(miRNA_pheno$new260) # More normal

hist(miRNA_pheno$RNA_quality)
miRNA_pheno$newquality <- ifelse(miRNA_pheno$RNA_quality < 4, NA, miRNA_pheno$RNA_quality)
hist(miRNA_pheno$newquality) # More normal

quantcon<-quantile(miRNA_pheno$concentration, probs = seq(0,1, 0.1))
miRNA_pheno$rankcon <- cut(miRNA_pheno$concentration, quantcon)
quantqual <-quantile(miRNA_pheno$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
miRNA_pheno$rankqual <- cut(miRNA_pheno$RNA_quality, quantqual)
quant260 <- quantile(miRNA_pheno$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
miRNA_pheno$rank260 <-cut(miRNA_pheno$`_260_280`, quant260)


#Spot Check miR-186-5p-a2
# Changes in MLR with added covariates
summary(lm(miR_186_5p_a2~concentration, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~`_260_280`, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~RNA_quality, data = miRNA_pheno))

summary(lm(miR_186_5p_a2~f8cbtobmd, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + concentration, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + RNA_quality, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + newquality, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankqual, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + `_260_280`, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + new260, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rank260, data = miRNA_pheno))
summary(lmer(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + (1|Isolation_Batch), data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con+ newquality, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con + newquality + `_260_280`, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con + newquality + new260, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual, data = miRNA_pheno))
summary(lm(miR_186_5p_a2~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno))


plot(miR_186_5p_a2~AGE8,data = miRNA_pheno)
abline(lm(miR_186_5p_a2~AGE8,data = miRNA_pheno), col = "red")
plot(miR_186_5p_a2~SEX,data = miRNA_pheno)
plot(miR_186_5p_a2~WGT8,data = miRNA_pheno)
abline(lm(miR_186_5p_a2~WGT8,data = miRNA_pheno), col = "red")
plot(miR_186_5p_a2~new260,data = miRNA_pheno)
abline(lm(miR_186_5p_a2~new260,data = miRNA_pheno), col = "red")
plot(miR_186_5p_a2~log_con,data = miRNA_pheno)
abline(lm(miR_186_5p_a2~log_con,data = miRNA_pheno), col = "red")


#Spot Check miR-19a-3p
# Changes in MLR with added covariates
summary(lm(miR_19a_3p~concentration, data = miRNA_pheno))
summary(lm(miR_19a_3p~`_260_280`, data = miRNA_pheno))
summary(lm(miR_19a_3p~RNA_quality, data = miRNA_pheno))

summary(lm(miR_19a_3p~f8cbtobmd, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + concentration, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + RNA_quality, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + newquality, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankqual, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + `_260_280`, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + new260, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rank260, data = miRNA_pheno))
summary(lmer(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + (1|Isolation_Batch), data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con+ newquality, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con + newquality + `_260_280`, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + log_con + newquality + new260, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual, data = miRNA_pheno))
summary(lm(miR_19a_3p~f8cbtobmd+AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno))

plot(miR_19a_3p~AGE8,data = miRNA_pheno)
abline(lm(miR_19a_3p~AGE8,data = miRNA_pheno), col = "red")
plot(miR_19a_3p~SEX,data = miRNA_pheno)
plot(miR_19a_3p~WGT8,data = miRNA_pheno)
abline(lm(miR_19a_3p~WGT8,data = miRNA_pheno), col = "red")
plot(miR_19a_3p~new260,data = miRNA_pheno)
abline(lm(miR_19a_3p~new260,data = miRNA_pheno), col = "red")
plot(miR_19a_3p~log_con,data = miRNA_pheno)
abline(lm(miR_19a_3p~log_con,data = miRNA_pheno), col = "red")

