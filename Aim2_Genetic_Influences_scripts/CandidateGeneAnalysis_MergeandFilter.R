# Read in all .Rdata files
rm(list = ls())
library(tidyverse)
library(xlsx)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")

genotype_files<- list.files("./", pattern = "_filtered_genotype.RData")
pheno <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/PhenoData_5_28.csv")

df <- NULL
for (i in genotype_files){
  print(i)
  load(i)
  gene <- merge(tidy$gt, tidy$fix, by = "POS")
  gene$Gene <- gsub(x = i, pattern = "_vcfr_tidy_filtered_genotype.RData",replacement ="")
  df <- rbind(df,gene)
}


# Filter for people with phenotype data and merge with pheno data
all_genes_pheno <- merge(df, pheno, by.x = "Indiv", by.y = "shareid")

# Number of SNPs (1482)
length(unique(all_genes_pheno$POS)) 

# Number of People (1527)
length(unique(all_genes_pheno$Indiv)) 

# Check that bad snps were filtered out
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
good_snps<-read.table("WellImputedPositions8.txt", header = T, sep = "\t")
bad_snps <- setdiff(unique(all_genes_pheno$POS),as.integer(good_snps$POS)) # Empty so no good snps made it this far

pheno_snp <- within(all_genes_pheno, Genotype <- relevel(Genotype, ref = "normal"))
pheno_snp_pos<- pheno_snp %>% group_by(.,POS,Gene) %>% nest()

table(pheno_snp_pos$Gene)

min(pheno_snp_pos$rows)
max(pheno_snp_pos$rows)
# No polyallelic SNPs present
save(pheno_snp_pos, file = "CandidateGenes_4_2020.RData")
