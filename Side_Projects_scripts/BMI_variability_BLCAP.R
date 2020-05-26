# BMI variability as a function of genotype of BLCAP SNPs
# SNPs: rs6090836 and rs6019102 and rs3764718
# Chromosome 20
# txt file with Chromosomes of interest
rm(list=ls())
library(dplyr)
library(vcfR)

options(stringsAsFactors = F)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")

# Write file of SNPs of Interest for VCF Extraction ####
rsIDs <- c("rs6090836", "rs6019102", "rs3764718")
pos <- c("36154004","36153344", "36138972")
SNPIds <- data.frame(rsID = rsIDs, POS = pos)
chrpos <- data.frame(CHR = 20, POS = pos)
write.table(chrpos,file = "SNPsOfInterest_BMIVariability.txt", quote = F, row.names = F, sep = '\t')

# Individuals of interest:####
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data/")
Cov <- read.delim("OffspringCovariates.txt", skip=10,header=T,stringsAsFactors = F)
write.table(Cov$shareid, "/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/BMI_Indiv.txt",
          quote = F, row.names = F, sep = '\t')

# Extract Covariates of Interest ####
pheno <- Cov %>% select(.,shareid, SEX, AGE8, BMI8, BG8, FASTING_BG8)
pheno1 <- na.omit(pheno)  
# Extracting SNPs portion from vcf files #### 
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")
system("vcftools --gzvcf chr20_c1.vcf.gz --chr 20 --positions SNPsOfInterest_BMIVariability.txt --keep BMI_Indiv.txt --recode --out BLCAPSNPsBMI")

# Read VCF File  to RData####
source("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Aim2_Genetic_Influences_scripts/vcftotidy_function.R")
vcftotidy("BLCAPSNPsBMI")
load("BLCAPSNPsBMI_vcfr_tidy_filtered_genotype.RData")

# Merge with Pheno Data and rsIDs
gene <- merge(tidy$gt, tidy$fix, by = "POS")
gene <- merge(gene, SNPIds, by = "POS")
all_SNPs_pheno <- merge(gene, pheno1, by.x = "Indiv", by.y = "shareid")
all_SNPs_pheno <- within(all_SNPs_pheno, Genotype <- relevel(Genotype, ref = "normal"))

# Select Columns of Interest ####
data <- all_SNPs_pheno %>% select(., Indiv, rsID, Genotype, BMI8, AGE8, SEX, BG8, FASTING_BG8)

# Convert from Long to Wide ####
all_SNPs_pheno_wide <- spread(data, key = rsID, value = Genotype)

write.csv(all_SNPs_pheno_wide, "/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data/BLCAP_BMI_SNPs.csv", quote = F, row.names = F)
