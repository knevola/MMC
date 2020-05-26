# SNP filtering
rm(list = ls())
library(dplyr)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
fem_results <- read.csv("AllSNPs_LMEKIN_results_female.csv")
male_results <- read.csv("AllSNPs_LMEKIN_results_male.csv")
fem_POS <-  c(115803375,239972561,240223080,106736732,60025809)
male_POS <- c(240050108,240112014,84682179,60001153,60026732,43177169)

fem_results_sig <- fem_results %>% filter(., POS %in% c(fem_POS, male_POS)) %>% filter(., var %in% c("gt_DS", "gt_DS:BBYes"))
male_results_sig <- male_results %>% filter(., POS %in% c(fem_POS,male_POS)) %>% filter(., var %in% c("gt_DS", "gt_DS:BBYes"))

                                            