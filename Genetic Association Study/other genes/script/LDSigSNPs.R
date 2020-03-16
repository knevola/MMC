# LD between Sig Snps
rm(list = ls())
library(gpart)
library(tidyverse)
library(genetics)
library(LDheatmap)
library(RColorBrewer)
library(gpart)
library(dplyr)
library(xlsx)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")
load("All_Genes_SNPs.RData")
all_genes <- saveLoadReference
Sig_snps <- read.csv("FrequencySigSNPs.csv")
gene_list <- read.xlsx("Candidate Gene List.xlsx", sheetIndex = 1)

genes_symbol<- as.character(gene_list$Gene.Symbol)
genes_symbol1 <- gsub("TNFS11", "RANKL", genes_symbol)
genes_symbol2 <- gsub("TNFRSF11A", "RANK", genes_symbol1)
genes_symbol2 <- gsub("TNFRSF11B", "OPG", genes_symbol2)
gene_list$Gene <- genes_symbol2

Sig_snps <- separate(Sig_snps, col = Var1, into = c("Gene", "Position"), sep = ":", remove = F)
snps_all <- merge(gene_list, Sig_snps, by = "Gene")

sig_Genes1<- all_genes %>% unnest(.,cols = c(data)) 
genotype <- sig_Genes1 %>% filter(., POS %in% snps_all$Position) %>% group_by(., POS) %>% nest()
geno1 <- NA
for (i in 1:length(genotype$POS)){
  print(i)
  g <- separate(genotype[[2]][[i]], col = gt_GT_alleles, into = c("a1", "a2"), sep = 2, remove = F)
  g$a3 <-substr(g$a1, start = 1, stop = 1)
  geno <-genotype(g$a3,g$a2)
  geno1 <-makeGenotypes(data.frame(geno1,geno))
}
geno1 <- geno1[-1]
    
LD_calc<- LD(g1 = geno1)
LD_pvalues<- LD_calc$`P-value`
names(LD_pvalues)<- genotype$POS
row.names(LD_pvalues)<- genotype$POS
LD_connected <- LD_pvalues < 0.05
LD_connected[is.na(LD_connected)] <- 1
amount_LD <- rowSums(LD_connected)
sum(amount_LD)
    
jpeg(paste("Pairwise LD for Sig SNPs.jpg"))
LDheatmap(geno1, genotype$POS, color = brewer.pal(n = 10,name = "RdBu"),SNP.name = genotype$POS, flip = T, title = "Pairwise LD for Sig SNPs")
dev.off()
write.csv(LD_pvalues, "PairwiseLDpvaluesSigSnps.csv", quote = F, row.names = T)
