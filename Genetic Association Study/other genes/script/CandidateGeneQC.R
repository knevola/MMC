rm(list = ls())
library(vcfR)
library(snpStats)
library(tidyverse)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes")
recode_vcf <- list.files(pattern = ".recode.vcf")
i = recode_vcf[1]
for (i in recode_vcf){
  system(paste("~/plink/plink --vcf ", i, " --biallelic-only list --out ", gsub(".recode.vcf", "", i), sep = ""))
}

# bed_files <- list.files(pattern = ".bed")
# bim_files <- list.files(pattern = ".bim")
# fam_files <- list.files(pattern = ".fam")
# Plink_files <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
# write.table(Plink_files, "Plink_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)
# 

system("~/plink/plink --merge-list Plink_files.txt --make-bed --out all_cand_genes_c1") # <- there are 23 variants that have 3 or more alleles present

genes <- gsub(pattern = ".bed", replacement = "", x = bed_files)

for (i in genes){
  system(paste("~/plink/plink --bfile ", i, " --exclude all_cand_genes_c1-merge.missnp --make-bed --out ", i, "-clean", sep = ""))
}

# Clean Files
bed_files <- list.files(pattern = "-clean.bed")
bim_files <- list.files(pattern = "-clean.bim")
fam_files <- list.files(pattern = "-clean.fam")
Plink_files <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
write.table(Plink_files, "Plink_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)

system("~/plink/plink --merge-list Plink_files.txt --make-bed --out all_cand_genes_c1")
system("~/plink/plink --bfile all_cand_genes_c1 --mind 0.10 --recode --out all_chr_c1_clean_mind") # Exclude samples missing 10% or more of genotype calls
system("~/plink/plink --file all_chr_c1_clean_mind --maf 0.05 --recode --out all_chr_c1_MAF_greater_5") # MAF > 0.05
system("~/plink/plink --file all_chr_c1_MAF_greater_5 --geno 0.05 --recode --out all_chr_c1_MAF_greater_5_clean") # Remove SNPs more than 5% missing genotypes

# Check duplicate samples/ relatedness (IBD)
system("~/plink/plink --file all_chr_c1_MAF_greater_5 --genome --out duplicates")

dups <- read.table("duplicates.genome", header = T)
problem_pairs = dups[which(dups$PI_HAT > 0.5),] # highly related individuals <- may want to repeat with larger number of SNPs

prob_people <- problem_pairs[1:2]
write.table(prob_people, "IBS_excluded.txt", sep = " ", row.names = F, quote = F, col.names = F)

system("~/plink/plink --file all_chr_c1_MAF_greater_5 --remove IBS_excluded.txt --make-bed --out all_cand_gene_clean3")
gene <- read.plink(bed = "all_cand_gene_clean3.bed",bim = "all_cand_gene_clean3.bim", fam = "all_cand_gene_clean3.fam" )
snp_list <- gene$map$position
id_list <- gene$fam$member

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
load("allSNPs_cleaned.RData")
length(intersect(tidy_gt_nest_clean$POS, snp_list)) # All SNPs filtered correctly

unnest_tidy <- tidy_gt_nest_clean %>% unnest(cols = c(data))
unnest_tidy <- within(unnest_tidy, Genotype <- relevel(Genotype, ref = "normal"))
tidy_clean2 <- unnest_tidy %>% filter(., Indiv %in% id_list) %>% group_by(., POS) %>% nest()
save(tidy_clean2, file = "allSNPs_cleaned_andQC.RData")
