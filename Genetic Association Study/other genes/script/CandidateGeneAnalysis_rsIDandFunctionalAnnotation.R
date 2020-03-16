rm(list = ls())
library(xlsx)
library(rentrez)
library(qqman)
library(tidyverse)
library(liftOver)
library(GEOquery)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

load("All_Genes_SNPs.RData")
all_snps <- saveLoadReference
gene_list <- read.xlsx("Candidate Gene List.xlsx", sheetIndex = 1)
all_snps_unnested <- all_snps %>% unnest(cols = c(data))
all_snps_names <- data.frame(Start = all_snps_unnested$POS, End = all_snps_unnested$POS, Gene= all_snps_unnested$Gene)
all_snps_names <- unique(all_snps_names)
genes_symbol<- as.character(gene_list$Gene.Symbol)
genes_symbol1 <- gsub("TNFS11", "RANKL", genes_symbol)
genes_symbol2 <- gsub("TNFRSF11A", "RANK", genes_symbol1)
genes_symbol2 <- gsub("TNFRSF11B", "OPG", genes_symbol2)
gene_list$Gene <- genes_symbol2
all_snps_names1 <- merge(all_snps_names, gene_list, by.x = "Gene", by.y = "Gene", all.x = T)
snp_names <- all_snps_names1[-c(4:12)]
snp_names$Start <- as.character(snp_names$Start)
snp_names$End <- as.character(snp_names$End)
snp_names$Chromosome <- as.character(snp_names$Chromosome)
chain<- import.chain("hg19ToHg38.over.chain")
SNPs<-makeGRangesFromDataFrame(snp_names)
seqlevelsStyle(SNPs) = "UCSC"  
SNPs38 = liftOver(SNPs, chain)


set_entrez_key("")
SNPs38_data <- SNPs38@unlistData
new_snps<- as.data.frame(SNPs38_data)

# Search Entrez for rsID
snps<-data.frame(SNPs = paste( gsub("chr", "",new_snps$seqnames),":",as.character(new_snps$start), sep = ""), Gene = snp_names$Gene, Old_Position = paste(snp_names$Chromosome, ":", snp_names$Start, sep = ""))

annotate_snp <- function(snps) {
  data = tibble(SNP = snps$SNPs, Gene = snps$Gene)
  data$rsID = NA
  data$effect = NA
  for (i in 1:length(data$SNP)) {
    starttime=Sys.time()  
    print(i)
    search=entrez_search("snp",snps$SNPs[i])
    Sys.sleep(0.3)
    search_id=entrez_summary(db="snp",id=search$ids)
    Sys.sleep(0.3)
    rs_id = NA
    fctn_class = NA
    if (length(search_id) == 31){
      rs_id = search_id$snp_id
      fctn_class = search_id$fxn_class
    } else {
      for (j in 1:length(search_id)){
        rs_id = rbind(rs_id,search_id[[j]][["snp_id"]])
        fctn_class = rbind(fctn_class, search_id[[j]][["fxn_class"]])
      }
      rs_id = rs_id[-1]
      fctn_class = fctn_class[-1]
    }
    data$rsID[i]=list(rs_id)
    data$effect[i]=list(fctn_class)
    Sys.sleep(3)
    print(Sys.time()-starttime)
    
  }
  return(data)
  
}
results <- annotate_snp(snps) # Will take about 4 to 5 hours to run.
results$old_POS <- paste(snp_names$Chromosome, ":", snp_names$Start, sep = "")
saveRDS(results, file = 'AllSNPsrsID.RData')
