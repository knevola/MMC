# QC on VCF Files
rm(list = ls())
library(GEOquery)
options(stringasFactor =F)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

# File with List of Well Imputed Snps
SNPs_file <- "WellImputedPositions8.txt"

# List of subjects
pheno <- read.csv("PhenoData_5_28.csv")
subjects <- pheno$shareid
write.table(subjects, "Subjects.txt", sep = '\t', row.names = F, col.names = F)
subject_file <- "Subjects.txt" 

# # Use vcf-concat to combine all vcf files and convert to PLINK OR
# system("vcf-concat *.vcf.gz > allchr_c1.vcf")
# # OR
# system("vcf-concat -f gzfile.txt > allchr_c1.vcf")
# 
# system("vcf-concat allchr_c1.vcf chr10_c1.vcf > allchr_c1.vcf")
# 
# system("~/plink/plink --vcf allchr_c1.vcf --out allchr_c1")

# Convert each VCF.gz to Plink and Merge Plink files
# gzfiles <- as.vector(as.matrix(read.table("gzfile.txt", sep = '\t', stringsAsFactors = F, header = F)))
# 
for (i in gzfiles){
  syst(paste("~/plink/plink --vcf ", i, "--extract ", SNPs_file," --make-bed --out", gsub(".vcf.gz", "", i), sep = ""))
}
system(paste('~/plink/plink --vcf chr10_c1.vcf ', "--extract ", subject_file,' --make-bed --out chr10_c1'))

#bed_files <- list.files("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data", pattern = ".bed")
#bim_files <- list.files("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data", pattern = ".bim")
#fam_files <- list.files("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data", pattern = ".fam")
#Plink_files <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
#write.table(bed_files, "bed_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#write.table(bim_files, "bim_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#write.table(fam_files, "fam_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#write.table(Plink_files, "Plink_files.txt", sep = "\t", row.names = F, col.names = F, quote = F)

system("~/plink/plink --merge-list Plink_files.txt --recode --out all_chr_c1")
                                                                                                                                                                                                                                                                                                                     #  How many variants are there in all_chr_c1

# Filter well imputed SNPs (at this point I should have a plink file called all_chr_c1)
system("~/plink/plink --file all_chr_c1.fam --mind 0.10 --recode --out all_chr_c1_clean_mind") # Exclude samples missing 10% or more of genotype calls
system(paste("~/plink/plink --file all_chr_c1_clean_mind --extract ", SNPs_file," --recode --out all_chr_c1_well_imputed", sep = "")) 

# Filter for subjects
system(paste("~/plink/plink --file all_chr_c1_well_imputed --extract ", subject_file, " --recode --out all_chr_c1_SNP_sub", sep = ""))

# Calculate and Filter MAF
system("~/plink/plink --file all_chr_c1_SNP_sub --maf 0.05 --recode --out all_chr_c1_MAF_greater_5") # MAF > 0.05
system("~/plink/plink --file all_chr_c1_MAF_greater_5 --geno 0.05 --recode --out all_chr_c1_MAF_greater_5_clean") # Remove SNPs more than 5% missing genotypes

# Check Sex and duplicate samples/ relatedness (IBD)
system("~/plink/plink --file all_chr_c1_MAF_greater 5 --check-sex --out all_chr_c1_sex_check.sexcheck")
sex_check <- read.table("all_chr_c1_sex_check.sexcheck", header = T)
sex_problem <- sex_check[sex_check$STATUS == "PROBLEM",]
sex_problem

system("~/plink/plink --file all_chr_c1_MAF_greater_5 --genome --out duplicates")
dups <- read.table("duplicates.genome", header = T)
problem_pairs <- dups[dups$PI_HAT > 0.4,]
hist(dups$PI_HAT)

#Create file of IBS_excluded individuals

# Filter out related individuals
system("~/plink/plink --file all_chr_c1_MAF_greater_5 --remove IBS_excluded.txt --recode --out all_chr_c1_clean2")
