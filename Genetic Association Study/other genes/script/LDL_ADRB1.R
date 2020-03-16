rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")
LDL_cov <- read.csv("LDL_cov.csv")
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
rsIDs<- readRDS("AllSNPsrsID.RData")
load("allSNPs_cleaned_andQC_genomewide.RData")

rsIDsADRB1 <- rsIDs[rsIDs$Gene == "ADRB1",]
rsIDsADRB1$POS <- gsub("10:","",rsIDsADRB1$old_POS)
ADRB1data <- tidy_clean2 %>% unnest(., cols = c(data)) %>% filter(.,Gene == "ADRB1")
rsID_merge <- merge(ADRB1data, rsIDsADRB1, by = "POS")

saveRDS(rsID_merge, "ADRB1_LDL_Data.RData")
readRDS("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ADRB1_LDL_Data.RData")
