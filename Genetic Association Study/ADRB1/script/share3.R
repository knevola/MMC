# Extract ADRB1 region from Chromosome 10 data
rm(list=ls())
library(vcfR)
library(R.utils)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/ADRB1/data")
#gunzip("chr10_c1.vcf.gz")
#get = read.vcfR(file="chr10_c1.vcf",cols=c(1:9))

# Region of ADRB1
min = 115803806 - 2000    # ADRB1 (uc001lba.3) at chr10:115803806-115806667 , 2kB upstream and 0.5 kB downstream
max = 115806667 + 500

# 
# # Hannahs Set numbers
# min = 115788936    # ADRB1 coords: 115,803,806 - 115,806,667
# max = 115821538

start = 1630256
stop = start + 439
adrb1 = read.vcfR(file="chr10_c1.vcf", skip = (start-1),nrows = (stop-start))

# ID numbers are not unique so set ID to arbitrary numbers
ids<-adrb1@fix[,3]
adrb1@fix[,3] <- 1:439
adrb1@fix[,3]

# Convert VcrF to tidy data
adrb1_tidy<- vcfR2tidy(adrb1, format_fields = c("GT", "DS", "GP"))

# Calculating the minor allele frequency
maf <- as.data.frame(maf(adrb1))
adrb1_tidy$meta
adrb1_tidy$fix
adrb1_tidy$gt
adrb1_tidy$fix$MAF <- maf$Frequency

save(adrb1_tidy, file = "adrb1_vcfr_tidy.RData")

# Pull out SNPs in ADBR1 gene or 2 kB upstream or 0.5 kB downstream of ADRB1 (60 SNPs)
adrb1_tidy$fix <- adrb1_tidy$fix[adrb1_tidy$fix$POS > min & adrb1_tidy$fix$POS < max,]
adrb1_tidy$gt <- adrb1_tidy$gt[adrb1_tidy$gt$POS > min & adrb1_tidy$gt$POS < max,]

# Filter for MAF > 0.05 (5 SNPs)
adrb1_tidy$fix <- adrb1_tidy$fix[adrb1_tidy$fix$MAF > 0.01,]
adrb1_tidy$gt <- adrb1_tidy$gt[adrb1_tidy$gt$POS %in% adrb1_tidy$fix$POS,]

# Pull out data for Individuals with PhenoData
pheno <- read.csv("PhenoData_5_28.csv")
adrb1_tidy$gt <- adrb1_tidy$gt[adrb1_tidy$gt$Indiv %in% pheno$shareid,]
length(unique(intersect(pheno$shareid,adrb1_tidy$gt$Indiv))) # 1527 of 1645 in pheno data have genotype data

save(adrb1_tidy, file = "adrb1_vcfr_tidy_filtered.RData")

