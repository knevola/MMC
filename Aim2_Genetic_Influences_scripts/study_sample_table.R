rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(gtsummary)

load("CandidateGenes_4_2020.RData")
pheno <- pheno_snp_pos$data[[1]]
pheno1 <- pheno[18:38]
pheno1 %>% tbl_summary(.,by = "SEX",statistic = list(all_continuous() ~ "{mean} ({sd})",
                                                     all_categorical() ~ "{n} / {N} ({p}%)")) 
pheno1 %>% tbl_summary(., statistic = list(all_continuous() ~ "{mean} ({sd})",
                                           all_categorical() ~ "{n} / {N} ({p}%)")) 
