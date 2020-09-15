rm(list=ls())
library(vcfR)
library(R.utils)
setwd("Z:/65488/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v18.p11.c1.HMB-IRB-MDS/GenotypeFiles/phg000835.v2.FHS_SHARE_imputed_HRC1.genotype-imputed-data.c1.HMB-IRB-MDS/phg000835.v2.FHS_SHARE_imputed_HRC1.genotype-imputed-data.c1")
#gunzip("chr5_c1.vcf.gz")
#get = read.vcfR(file="chr5_c1.vcf",cols=c(1:9))
min = 148186349
max = 148828687
hg38start = 148825245
start = 2109948
stop = start + 9486
adrb2 = read.vcfR(file="chr5_c1.vcf",skip = (start-1),nrows = (stop-start))
samples = dimnames(adrb2@gt)[[2]]
cohort = read.delim(file="cohort.txt",header=T,stringsAsFactors=F)
shareid_int = as.integer(samples[2:length(samples)])
length(which(shareid_int %in% cohort$shareid)) # 1925 out of 2071
pull = adrb2@fix[,8]
ac = strsplit(pull,split=";")
num_alt = array(0,dim=c(length(ac)))
for (i in 1:length(ac)) {
  pull2 = strsplit(ac[[i]][1],split="=")
  num_alt[i] = as.numeric(pull2[[1]][2])
}
maf = num_alt/14324
hist(maf)
keep = which(maf>0.05)
adrb2_reduced = adrb2
adrb2_reduced@fix = adrb2@fix[keep,]
adrb2_reduced@gt = adrb2@gt[keep,]
write.vcf(adrb2_reduced,"adrb2.vcf")
check = read.vcfR(file="adrb2.vcf")
