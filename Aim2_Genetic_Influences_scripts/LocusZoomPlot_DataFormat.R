# LocusZoom Plot Data Format
# Can upload chromosome, position, reference and alt alleles, and p-value for entire analysis to web to make interactive plots
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(dplyr)
library(tidyr)

# Load Data
fem_results <- read.csv("AllSNPs_LMEKIN_results_female.csv")
male_results <- read.csv("AllSNPs_LMEKIN_results_male.csv")
load('AllSNPsrsID_4_2020.RData')
load("CandidateGenes_4_2020.RData")

# Filter for data of interest
fem_result_int <- fem_results %>% filter(., var == "gt_DS:BBYes")
male_results_int <- male_results %>% filter(., var == "gt_DS:BBYes")

# Separate Position marker
results <- results %>% separate(., old_POS, into = c("Old_Chr", "Old_Pos"), sep = ":", remove = F)
results <- results %>% separate(., SNP, into = c("Chromosome", "New_Pos"), sep = ":", remove = F)

# Merge results with positions
fem_result_int_info <- merge(fem_result_int, results, by.x = "POS", by.y = "Old_Pos")
male_result_int_info <- merge(male_results_int, results, by.x = "POS", by.y = "Old_Pos")

# Extract reference and alterative alleles
data <- pheno_snp_pos %>% unnest(cols = c(data))
data_info <- data %>% group_by(., POS, gt_GT_alleles, Genotype, MAF) %>% nest()
data_info <- data_info[1:4]

data_info1 <- data_info[data_info$Genotype %in% c("normal", "homozygous alternative"),]
data_info1$allele <- substr(data_info1$gt_GT_alleles, start = 1, stop = 1)

data_info2 <- data_info1[-2]
data_info2 <- spread(data_info2, Genotype, allele)
names(data_info2)<- c("POS", "MAF","ref_allele", "alt_allele")

# Merge alleles into results
fem_result_int_info1 <- merge(fem_result_int_info, data_info2, by = "POS")
fem_result_int_info1$Old_POS_Chr <- paste(fem_result_int_info1$Chromosome, fem_result_int_info1$POS, sep = ":")
male_result_int_info1 <- merge(male_result_int_info, data_info2, by = "POS")
male_result_int_info1$Old_POS_Chr <- paste(male_result_int_info1$Chromosome, male_result_int_info1$POS, sep = ":")

# Select columns of interest
fem_LZP_data <- fem_result_int_info1 %>% dplyr::select(Chromosome, New_Pos, ref_allele, alt_allele, p, Gene.y)
male_LZP_data <- male_result_int_info1 %>% dplyr::select(Chromosome, New_Pos, ref_allele, alt_allele, p, Gene.y)

write.table(fem_LZP_data, "Fem_LocusZoomPlot_data.txt", quote = F, sep = "\t", row.names = F)
write.table(male_LZP_data, "Male_LocusZoomPlot_data.txt", quote = F, sep = "\t", row.names = F)

# Select data for GCTA COJO
fem_GCTA <- fem_result_int_info1 %>% dplyr::select(Old_POS_Chr,alt_allele, ref_allele,MAF, Value, Std.Error,p,n )
male_GCTA <- male_result_int_info1 %>% dplyr::select(Old_POS_Chr, alt_allele, ref_allele, MAF, Value, Std.Error, p, n)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")
write.table(fem_GCTA, "GCTA_fem.ma", sep = " ", row.names = F, quote = F)
write.table(male_GCTA, "GCTA_male.ma", sep = " ", row.names = F, quote = F)
