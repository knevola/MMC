# Query Biomart
rm(list = ls())
library(xlsx)
library(tidyverse)
library(biomaRt)
library(liftOver)
library(rentrez)

Sig_snps <- read.csv("FrequencySigSNPs.csv")
gene_list <- read.xlsx("Candidate Gene List.xlsx", sheetIndex = 1)

genes_symbol<- as.character(gene_list$Gene.Symbol)
genes_symbol1 <- gsub("TNFS11", "RANKL", genes_symbol)
genes_symbol2 <- gsub("TNFRSF11A", "RANK", genes_symbol1)
genes_symbol2 <- gsub("TNFRSF11B", "OPG", genes_symbol2)
gene_list$Gene <- genes_symbol2

Sig_snps <- separate(Sig_snps, col = Var1, into = c("Gene", "Position"), sep = ":", remove = F)
snps_all <- merge(gene_list, Sig_snps, by = "Gene")

# Lift over to new genome build
snps_BM <- data.frame(CHR = snps_all$Chromosome, POS = snps_all$Position)
chain<- import.chain("hg19ToHg38.over.chain")
SNPs<-makeGRangesFromDataFrame(data.frame(snps_BM, Start = as.character(snps_BM$POS), End = as.character(snps_BM$POS)))
seqlevelsStyle(SNPs) = "UCSC"  
SNPs38 = liftOver(SNPs, chain)
SNPs38_data <- SNPs38@unlistData
new_snps<- as.data.frame(SNPs38_data)

# Search Entrez for rsID
snps<-data.frame(SNPs = paste( gsub("chr", "",new_snps$seqnames),":",as.character(new_snps$start), sep = ""), Gene = snps_all$Gene)

annotate_snp <- function(snps) {
  data = data.frame(SNP = snps$SNPs, Gene = snps$Gene)
  data$rsID = NA
  data$pubmed  = NA
  data$effect = NA  
  
  for (i in 1:length(data$SNP)) {
    starttime=Sys.time()  
    print(i)
    search=entrez_search("snp",snps$SNPs[i])
    Sys.sleep(0.3)
    search_id=entrez_summary(db="snp",id=search[["ids"]][[1]])
    Sys.sleep(0.3)
    rs_id = search_id[["snp_id"]]
    
    data$rsID[i]=paste("rs",rs_id, sep = '')
    
    pubmed_id=entrez_link(dbfrom="snp",db="pubmed",id=rs_id)
    
    
    link = ""
    Sys.sleep(0.3)
    if (length(pubmed_id[["links"]][["snp_pubmed"]])>0) {
      for (j in 1:length(pubmed_id[["links"]][["snp_pubmed"]])) {
        
        pubmedlink = entrez_link(dbfrom="pubmed",id=as.numeric(pubmed_id[["links"]][["snp_pubmed"]][[j]]),cmd="llinks")
        Sys.sleep(0.3)
        link=paste(link,pubmedlink[["linkouts"]][[1]][[1]][["Url"]],sep = ",")
      }
      data$pubmed[i]=link
    }
    print(Sys.time()-starttime)
    
    
    
    data$effect[[i]]=search_id$fxn_class
    
  }
  
  return(data)
  
}
results <- annotate_snp(snps)
write.xlsx(results, "rsIDandAnnotationofSigSNPs.xlsx", row.names = F)
