# JCEM Meta-Analysis
rm(list = ls())
setwd("~/Thesis")
library(tidyverse)
library(xlsx)
library(meta)
meta_analysis<-read.xlsx("GeneticMetaAnalysis.xlsx", sheetIndex = 1)

ADRB1 <- meta_analysis %>% filter(., Gene == "ADRB1")
HDAC4 <- meta_analysis %>% filter(., Gene == "HDAC4")
RANK1 <- meta_analysis %>% filter(., rsID == "rs34170507")
RANK2 <- meta_analysis %>% filter(., rsID == "rs6567268")

ADRB1_meta <- metagen(TE = ADRB1$beta, seTE = ADRB1$se, studlab = ADRB1$Cohort,pval = ADRB1$p)
summary(ADRB1_meta)
forest(ADRB1_meta, digits = 4)

HDAC4_meta <- metagen(TE = HDAC4$beta, seTE = HDAC4$se, studlab = HDAC4$Cohort,pval = HDAC4$p)
summary(HDAC4_meta)
forest(HDAC4_meta, digits = 4)
