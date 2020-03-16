##Identify clusters of SNPs in LD for each gene:
rm(list = ls())
library("VariantAnnotation")
library("adjclust")
library(GGtools)
library("gpart")

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes")
files <-list.files(pattern = "genotype.RData")

#Gene list: 
candidate = gsub(pattern = "_.*?$", replacement = "",x = files)
output = data.frame(candidate, total.hypotheses = rep(0,length(candidate)))
#Iterate through the list
for (i in 1:length(candidate)) {
  print(candidate[i])
  #First input a genotype matrix:
  load(files[i])
  #Create a SNPinfo dataframe:
  SNPinfo = data.frame(tidy$fix$CHROM,as.character(tidy$fix$POS),as.numeric((tidy$fix$POS)))
  colnames(SNPinfo)=c("chrN","rsID","bp")
  #Create a genotype matrix : 
  geno = matrix(tidy$gt$gt_GT, ncol =length(tidy$fix$CHROM),byrow =T)
  geno[geno=="0|1" | geno =="1|0"]=1
  geno[geno=="0|0"]=2
  geno[geno=="1|1"]=0
  
  mode(geno)="numeric"
  geno=as.data.frame(geno)
  colnames(geno)=unique(tidy$gt$POS)
  #Pull out the list of bins size for each SNP cluster, based on D': 
  clq = CLQD(geno,SNPinfo,LD="Dprime")
  sum_hypotheses = length(table(clq))+sum(is.na(clq))
  print(paste("For", candidate[i],"the total is", as.integer(sum_hypotheses)))
  output$total.hypotheses[[i]]=sum_hypotheses
}
write.csv(output,file="/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/cluster.csv")
