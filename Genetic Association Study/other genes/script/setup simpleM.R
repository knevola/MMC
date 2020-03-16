#Getting SNP Sample Data Ready for simpleM
#example SNP file format:
  # row => SNPs
  # column => Unrelated individuals 
  
  # The data file should contain only POLYMORPHIC SNPs. 
  
  # Missing values should be imputed. 
  # There should be NO missing values in the SNP data file.
  # SNPs are coded as 0, 1 and 2 for the number of reference alleles. 
  # SNPs are separated by one-character spaces. 
  
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")
library(dplyr)
library(tidyverse)
library(reshape2)

`%notin%` <- Negate(`%in%`)
load("All_Genes_SNPs.RData")
all_snps <- saveLoadReference
tidy_gt <- all_snps %>% unnest(cols = c(data))
tidy_gt$simpleM[tidy_gt$Genotype == "normal"] <- 2 
tidy_gt$simpleM[tidy_gt$Genotype == "homozygous alternative"] <- 0 
tidy_gt$simpleM[tidy_gt$Genotype == "heterozygous"] <- 1 
tidy_gt$Indiv <- as.factor(tidy_gt$Indiv)
tidy_gt$POS <- as.factor(tidy_gt$POS)

wide <-dcast(tidy_gt, POS ~ Indiv, value.var = "simpleM", mean)
wide1 <-dcast(tidy_gt, POS ~ Indiv, value.var = "simpleM", length)
df3 <- data.frame(POS= wide1$POS, sum = rowSums(wide1[-1]))
bad_snps <- df3 %>% filter(., sum > 1527)
bad1 <- as.integer(as.character(bad_snps$POS))
tidy_gt$POS <- as.integer(as.character(tidy_gt$POS))
tidy_gt_clean <- tidy_gt %>% filter(., POS %notin% bad1)
write.csv(wide1, "duplicatesamplescheck.csv",quote = F, row.names = F)
wide <-dcast(tidy_gt_clean, POS ~ Indiv, value.var = "simpleM", mean)
tidy_gt_nest_clean<- tidy_gt_clean %>% group_by(., POS) %>% nest()
save(tidy_gt_nest_clean, file = "allSNPs_cleaned.RData")
# Issue: some individuals have multiple calls for snps, so used mean, still have issue if 1 and 2 or 0 and 2
wide_data <- wide[-1]
# ceiling_wide <- ceiling(wide_data)
# floor_wide <- floor(wide_data)

write.table(wide_data, file = "simpleM_allSNPs.txt", sep = " ", quote = F, row.names = F, col.names = F)
# write.table(ceiling_wide, file = "simpleM_allSNPs_ceiling.txt", sep = " ", quote = F, row.names = F, col.names = F)
# write.table(floor_wide, file = "simpleM_allSNPs_floor.txt", sep = " ", quote = F, row.names = F, col.names = F)
