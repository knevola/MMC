# SEX Stratified Model for top 2 miRNA
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

pheno_F <- miRNA_pheno[miRNA_pheno$SEX == "Female",]
pheno_M <- miRNA_pheno[miRNA_pheno$SEX == "Male",]

Fto19F<-summary(lm(miR_19a_3p ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
Fto19F<-as.data.frame(Fto19F$coefficients)
Fto19F$Group <- "Female"
Fto19M<-summary(lm(miR_19a_3p ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
Fto19M<-as.data.frame(Fto19M$coefficients)
Fto19M$Group <- "Male"
Fto19<-summary(lm(miR_19a_3p ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
Fto19<-as.data.frame(Fto19$coefficients)
Fto19$Group <- "Both"
LS19F<-summary(lm(miR_19a_3p ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
LS19F <- as.data.frame(LS19F$coefficients)
LS19F$Group <- "Female"
LS19M <- summary(lm(miR_19a_3p ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
LS19M <- as.data.frame(LS19M$coefficients)
LS19M$Group <- "Male"
LS19<-summary(lm(miR_19a_3p ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
LS19<-as.data.frame(LS19$coefficients)
LS19$Group <- "Both"
BB19F <- summary(lm(miR_19a_3p ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
BB19F <- as.data.frame(BB19F$coefficients)
BB19F$Group <- "Female"
BB19M <- summary(lm(miR_19a_3p ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
BB19M <- as.data.frame(BB19M$coefficients)
BB19M$Group <- "Male"
BB19 <- summary(lm(miR_19a_3p ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
BB19 <- as.data.frame(BB19$coefficients)
BB19$Group <- "Both"

Fto186F<-summary(lm(miR_186_5p_a2 ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
Fto186F<-as.data.frame(Fto186F$coefficients)
Fto186F$Group <- "Female"
Fto186M<-summary(lm(miR_186_5p_a2 ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
Fto186M<-as.data.frame(Fto186M$coefficients)
Fto186M$Group <- "Male"
Fto186<-summary(lm(miR_186_5p_a2 ~ f8cbtobmd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
Fto186<-as.data.frame(Fto186$coefficients)
Fto186$Group <- "Both"
LS186F<-summary(lm(miR_186_5p_a2 ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
LS186F <- as.data.frame(LS186F$coefficients)
LS186F$Group <- "Female"
LS186M <- summary(lm(miR_186_5p_a2 ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
LS186M <- as.data.frame(LS186M$coefficients)
LS186M$Group <- "Male"
LS186<-summary(lm(miR_186_5p_a2 ~ s8cbl24bd + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
LS186<-as.data.frame(LS186$coefficients)
LS186$Group <- "Both"
BB186F <- summary(lm(miR_186_5p_a2 ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_F))
BB186F <- as.data.frame(BB186F$coefficients)
BB186F$Group <- "Female"
BB186M <- summary(lm(miR_186_5p_a2 ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = pheno_M))
BB186M <- as.data.frame(BB186M$coefficients)
BB186M$Group <- "Male"
BB186 <- summary(lm(miR_186_5p_a2 ~ BB + AGE8 + HGT8 + WGT8 + rankcon + rank260 + rankqual,data = miRNA_pheno))
BB186 <- as.data.frame(BB186$coefficients)
BB186$Group <- "Both"

MLR19<-rbind(Fto19[2,], Fto19F[2,], Fto19M[2,], LS19[2,], LS19F[2,], LS19M[2,], BB19[2,], BB19F[2,], BB19M[2,])
MLR186<-rbind(Fto186[2,], Fto186F[2,], Fto186M[2,], LS186[2,], LS186F[2,], LS186M[2,], BB186[2,], BB186F[2,], BB186M[2,])

write.csv(MLR19, "SexStratifiedMiR19.csv")
write.csv(MLR186, "SexStratifiedMiR186.csv")
