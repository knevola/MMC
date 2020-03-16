rm(list = ls())
library(vcfR)
library(snpStats)
library(tidyverse)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")


bed_files <- list.files(pattern = "_c1.bed")
bim_files <- list.files(pattern = "_c1.bim")
fam_files <- list.files(pattern = "_c1.fam")
fam_files <- fam_files[-1]
Plink_files <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
write.table(Plink_files, "Plink_files_c1.txt", sep = "\t", row.names = F, col.names = F, quote = F)
# 

system("~/plink/plink --merge-list Plink_files_c1.txt --make-bed --out all_chr_c1") 
genes <- gsub(pattern = ".bed", replacement = "", x = bed_files)

for (i in genes){
  system(paste("~/plink/plink --bfile ", i, " --exclude all_chr_c1-merge.missnp --keep-fam Subjects.txt --make-bed --out ", i, "-clean", sep = ""))
}


# Clean Files
bed_files <- list.files(pattern = "-clean.bed")
bim_files <- list.files(pattern = "-clean.bim")
fam_files <- list.files(pattern = "-clean.fam")
Plink_files <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
write.table(Plink_files, "Plink_files_clean.txt", sep = "\t", row.names = F, col.names = F, quote = F)

system("~/plink/plink --merge-list Plink_files_clean.txt  --make-bed --out all_chr_c1_clean")
system("~/plink/plink --bfile all_chr_c1_clean --mind 0.10 --recode --out all_chr_c1_clean_mind") # Exclude samples missing 10% or more of genotype calls

system("~/plink/plink --file all_chr_c1_clean_mind --maf 0.05 --recode --out all_chr_c1_MAF_greater_5") # MAF > 0.05
system("~/plink/plink --file all_chr_c1_MAF_greater_5 --geno 0.05 --recode --out all_chr_c1_MAF_greater_5_clean") # Remove SNPs more than 5% missing genotypes

# Check duplicate samples/ relatedness (IBD)
system("~/plink/plink --file all_chr_c1_MAF_greater_5_clean --genome --out duplicates")

rm(list = ls())
dups <- read.table("duplicates.genome", header = T)
problem_pairs = dups[which(dups$PI_HAT > 0.4),] # highly related individuals <- may want to repeat with larger number of SNPs

prob_people <- problem_pairs[1:2]
write.table(prob_people, "IBS_excluded_genomewide.txt", sep = " ", row.names = F, quote = F, col.names = F)

system("~/plink/plink --file all_chr_c1_MAF_greater_5 --remove IBS_excluded_genomewide.txt --make-bed --out all_chr_clean3")
gene <- read.plink(bed = "all_chr_clean3.bed",bim = "all_chr_clean3.bim", fam = "all_chr_clean3.fam" )
id_list <- gene$fam$member

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
load("allSNPs_cleaned.RData")

unnest_tidy <- tidy_gt_nest_clean %>% unnest(cols = c(data))
unnest_tidy <- within(unnest_tidy, Genotype <- relevel(Genotype, ref = "normal"))
tidy_clean2 <- unnest_tidy %>% filter(., Indiv %in% id_list) %>% group_by(., POS) %>% nest()
save(tidy_clean2, file = "allSNPs_cleaned_andQC_genomewide.RData")
