rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")

# Convert to Plink
recode_vcf <- list.files(pattern = ".recode.vcf")
for (i in recode_vcf){
  system(paste("~/plink/plink --vcf ", i, " --biallelic-only list --out ", gsub(".recode.vcf", "", i), sep = ""))
}

genes <- gsub(pattern = ".recode.vcf", replacement = "", x = recode_vcf)
genes <- genes[-c(3,8)]

for (i in genes){
  system(paste("~/gcta_1.93.0beta/gcta64 --bfile ", i, " --cojo-file GCTA_fem.ma --cojo-slct --cojo-p 0.05 --diff-freq 1 --out ", i, "fem", sep = "" ))
}

for (i in genes){
  system(paste("~/gcta_1.93.0beta/gcta64 --bfile ", i, " --cojo-file GCTA_male.ma --cojo-slct --cojo-p 0.05 --diff-freq 1 --out ", i, "male", sep = "" ))
}

fem_results <- list.files(pattern = "fem.jma.cojo")
male_results <- list.files(pattern = "male.jma.cojo")

results <- c(fem_results, male_results)
data1 <- NULL

for (i in results){
  data <- read.table(i, header = T)
  data$file <- gsub(pattern = ".jma.cojo", replacement = "",x = i)
  data1 <- rbind(data1, data)
}

write.csv(data1, "Significant_SNPs_GCTA_COJO.csv", row.names = F, quote = F)
