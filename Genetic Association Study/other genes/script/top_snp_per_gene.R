rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
all_mlr <- read.csv("AllMLR_Results.csv", header = T)
all_mlr_LSFN <- all_mlr %>%  filter(., Model %in% c("FN:F", "FN:SC", "LS:F", "LS:SC")) %>% 
  filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
gene4_mlr <- all_mlr_LSFN %>% filter(., Pr...t.. < 0.0008333333)
all_mlr_LSFN_genes <- all_mlr_LSFN %>% group_by(Gene) %>% top_n(n = -1, wt = Pr...t..)
table(all_mlr_LSFN_genes$Gene)
write.csv(all_mlr_LSFN_genes, "TopPerSNP1.csv")

top_snps<- read.csv("TopPerSNP.csv")
snps_fa<- readRDS("AllSNPsrsID.RData")
snps_fa$effect <- as.character(snps_fa$effect)
snps_fa$rsID <- as.character(snps_fa$rsID)
top_snps1 <- top_snps[-c(6,15:23,25:36),]
top_snps_fa <- snps_fa %>% filter(., old_POS %in% top_snps1$SNP)
top_snps_fa$rsID <- as.character(top_snps_fa$rsID)
top_snps_fa$rsID <- gsub("c\\(", "", top_snps_fa$rsID)
top_snps_fa$rsID <- gsub("\\)", "", top_snps_fa$rsID)

top_snps_fa$effect <- as.character(top_snps_fa$effect)
top_snps_fa$effect <- gsub("c\\(", "", top_snps_fa$effect)
top_snps_fa$effect <- gsub("\\)", "", top_snps_fa$effect)
top_snps_fa$effect <- gsub("_", " ", top_snps_fa$effect)
top_snps_fa$effect <- noquote(top_snps_fa$effect)

write.table(top_snps_fa, "TopSNPsFA.txt", quote = F, row.names = F, sep = "\t")

table(top_snps1$Model)
table(top_snps1$var,top_snps1$Model)
snps_FA_mis <- snps_fa[grepl(pattern = "missense",x = snps_fa$effect),]

missense<- snps_FA_mis %>% separate(., col = old_POS, into = c("Chr","Position"), sep = ":",remove = F)
missense$POS <- paste(missense$Gene, missense$Position, sep = ":")

missense_models <- all_mlr_LSFN %>% filter(., SNP %in% missense$POS)

missense_gene4 <- missense_models %>% filter(., Pr...t.. < 0.0008333333) 
missense_gene <- missense_models %>% filter(., Pr...t.. < 0.003333333) 
sig_missense<- merge(missense_gene, missense, by.x = "SNP", by.y = "POS")
sig_missense$rsID <- as.character(sig_missense$rsID)
sig_missense$rsID <- gsub("c\\(", "", sig_missense$rsID)
sig_missense$rsID <- gsub("\\)", "", sig_missense$rsID)

sig_missense$effect <- as.character(sig_missense$effect)
sig_missense$effect <- gsub("c\\(", "", sig_missense$effect)
sig_missense$effect <- gsub("\\)", "", sig_missense$effect)
sig_missense$effect <- gsub("_", " ", sig_missense$effect)
sig_missense$effect <- noquote(sig_missense$effect)
write.table(sig_missense, "SignificantMissense.txt", sep = '\t', quote = F, row.names = F)

write.csv(top_snps1, "TopSNPsfinallist.csv", quote = F, row.names = F)
top_snps1$SNP<- gsub("^..:",":",top_snps1$SNP) 
top_snps1$SNP <- paste(top_snps1$Gene, top_snps1$SNP, sep="")
validSNPs <- all_mlr_LSFN %>% filter(., SNP %in% top_snps1$SNP)
write.csv(validSNPs, "SNPsforValidatonResults.csv", quote = F, row.names = F)
