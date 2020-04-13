# extracting gene portion from vcf files: 
#Read gene list:
rm(list=ls())
library(dplyr)
library(xlsx)
options(stringsAsFactors = F)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")
gene_list1 <- read.xlsx("CandidateGeneList_4_10.xlsx", sheetIndex = 1)
gene_list <- gene_list1 %>% filter(., Selected.for.analysis == "X")

#####extract gene from the chromosome: 
extract <-function(gene_list) {
  for (i in 1:nrow(gene_list)) {
    if (gene_list$Chromosome == 10){
      print(i)
      cmd = "vcftools --vcf chr"
      cmd = paste(cmd,gene_list$Chromosome[i],"_c1.vcf", " --chr ",gene_list$Chromosome[i]," --from-bp ",gene_list$Min[i]," --to-bp ",gene_list$Max[i], sep = "")
      cmd = paste(cmd," --positions WellImputedPositions8.txt --keep IndivWithPheno.txt --recode ", "--out ",gene_list$Gene.Symbol[i],sep="")
      system(cmd)
    }else{
      print(i)
      cmd = "vcftools --gzvcf chr"
      cmd = paste(cmd,gene_list$Chromosome[i],"_c1.vcf.gz", " --chr ",gene_list$Chromosome[i]," --from-bp ",gene_list$Min[i]," --to-bp ",gene_list$Max[i], sep = "")
      cmd = paste(cmd," --positions WellImputedPositions8.txt --keep IndivWithPheno.txt --recode ", "--out ",gene_list$Gene.Symbol[i],sep="")
      system(cmd)
    }
  }
}

extract(gene_list = gene_list) # This will take about 3.5 hours to run

print("Done.")

#chr 10 was loaded as vcf

i = 12
cmd = "vcftools --vcf chr"
cmd = paste(cmd,gene_list$Chromosome[i],"_c1.vcf", " --chr ",gene_list$Chromosome[i]," --from-bp ",gene_list$Min[i]," --to-bp ",gene_list$Max[i], sep = "")
cmd = paste(cmd," --positions WellImputedPositions8.txt --keep IndivWithPheno.txt --recode ", "--out ",gene_list$Gene.Symbol[i],sep="")
system(cmd)
