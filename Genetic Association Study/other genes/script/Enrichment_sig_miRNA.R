rm(list = ls())
library(tidyverse)
library(multiMiR)
library(edgeR)
setwd("/home/clary@mmcf.mehealth.org/Framingham/eQTL Analysis/data")
miRNAs <- read.csv("SignificantMiRNAswithGenes.csv", stringsAsFactors = F)
miRNAs$miRNA <- as.character(miRNAs$Var2)
miRNAs$miRNA <- gsub(x = miRNAs$miRNA,pattern ="_",replacement = "-" )
miRNAs$miRNA <- gsub(x = miRNAs$miRNA, pattern = "miR", replacement = "hsa-miR")
miRNAs$miRNA <- gsub(x = miRNAs$miRNA, pattern = "5p-a2", replacement = "5p")
miRNAs$miRNA <- gsub(x = miRNAs$miRNA, pattern = "5p-a1", replacement = "5p")

# miRNA per SNP
RANK_miRNA <- miRNAs[miRNAs$Var1 == "RANK",]
ADRB2_miRNA <- miRNAs[miRNAs$Var1 == "ADRB2",]
FoxO1_miRNA <- miRNAs[miRNAs$Var1 == "FoxO1",]

# Gene Targets of miRNA
genes<- get_multimir(mirna = unique(miRNAs$miRNA), summary = TRUE)
RANK_genes <- get_multimir(mirna = unique(RANK_miRNA$miRNA), summary = TRUE)
ADRB2_genes <- get_multimir(mirna = unique(ADRB2_miRNA$miRNA), summary = TRUE)
FoxO1_genes <- get_multimir(mirna = unique(FoxO1_miRNA$miRNA), summary = TRUE)

genesdf <- genes@data
RANK_genesdf <- RANK_genes@data
ADRB2_genesdf <- ADRB2_genes@data
FoxO1_genesdf <- FoxO1_genes@data

# Which miRNA are not found
setdiff(miRNAs$miRNA,unique(genesdf$mature_mirna_id))

genesdfu <- genesdf[3:5]
RANK_genesdfu <- RANK_genesdf[3:5]
ADRB2_genesdfu <- ADRB2_genesdf[3:5]
FoxO1_genesdfu <- FoxO1_genesdf[3:5]

genesdfuni <- unique(genesdfu)
RANK_genesdfuni <- unique(RANK_genesdfu)
ADRB2_genesdfuni <- unique(ADRB2_genesdfu)
FoxO1_genesdfuni <- unique(FoxO1_genesdfu)

length(unique(genesdf$target_entrez))
length(unique(RANK_genesdf$target_entrez))
length(unique(ADRB2_genesdf$target_entrez))
length(unique(FoxO1_genesdf$target_entrez))

overlapping_targets<- intersect(intersect(RANK_genesdf$target_symbol, ADRB2_genesdf$target_symbol), FoxO1_genesdf$target_symbol)

table<-as.data.frame(table(genesdfuni$target_entrez))
RANK_table <- as.data.frame(table(RANK_genesdfuni$target_entrez)) 
ADRB2_table <- as.data.frame(table(ADRB2_genesdfuni$target_entrez)) 
FoxO1_table <- as.data.frame(table(FoxO1_genesdfuni$target_entrez)) 

table2<-table[table$Freq >= 2 & table$Freq < max(table$Freq),]
RANK_table2<-RANK_table[RANK_table$Freq >= 2 & RANK_table$Freq < max(RANK_table$Freq),]
ADRB2_table2<-ADRB2_table[ADRB2_table$Freq >= 2 & ADRB2_table$Freq < max(ADRB2_table$Freq),]
FoxO1_table2<-FoxO1_table[FoxO1_table$Freq >= 2 & FoxO1_table$Freq < max(FoxO1_table$Freq),]

overlap_table2 <- intersect(RANK_table2$Var1, intersect(ADRB2_table2$Var1, FoxO1_table2$Var1))

table3<-table[table$Freq >= 3 & table$Freq < max(table$Freq),]
RANK_table3<-RANK_table[RANK_table$Freq >= 3 & RANK_table$Freq < max(RANK_table$Freq),]
ADRB2_table3<-ADRB2_table[ADRB2_table$Freq >= 3 & ADRB2_table$Freq < max(ADRB2_table$Freq),]
FoxO1_table3<-FoxO1_table[FoxO1_table$Freq >= 3 & FoxO1_table$Freq < max(FoxO1_table$Freq),]

overlap_table3 <- intersect(RANK_table3$Var1, intersect(ADRB2_table3$Var1, FoxO1_table3$Var1))

overlap<- genesdfuni %>% filter(., target_entrez %in% overlap_table3)
length(unique(overlap$mature_mirna_id))


GO3 <- goana(as.vector(table3$Var1))
RANK_GO3 <- goana(as.vector(RANK_table3$Var1))
ADRB2_GO3 <- goana(as.vector(ADRB2_table3$Var1))
FoxO1_GO3 <- goana(as.vector(FoxO1_table3$Var1))

KEGG3 <- kegga(as.vector(table3$Var1))
RANK_KEGG3 <- kegga(as.vector(RANK_table3$Var1))
ADRB2_KEGG3 <- kegga(as.vector(ADRB2_table3$Var1))
FoxO1_KEGG3 <- kegga(as.vector(FoxO1_table3$Var1))

GO3_p <- GO3[GO3$P.DE < 0.05,]
RANK_GO3_p <- RANK_GO3[RANK_GO3$P.DE < 0.05,]
ADRB2_GO3_p <- ADRB2_GO3[ADRB2_GO3$P.DE < 0.05,]
FoxO1_GO3_p <- FoxO1_GO3[FoxO1_GO3$P.DE < 0.05,]

KEGG3_p <- KEGG3[KEGG3$P.DE < 0.05,]
RANK_KEGG3_p <- RANK_KEGG3[RANK_KEGG3$P.DE < 0.05,]
ADRB2_KEGG3_p <- ADRB2_KEGG3[ADRB2_KEGG3$P.DE < 0.05,]
FoxO1_KEGG3_p <- FoxO1_KEGG3[FoxO1_KEGG3$P.DE < 0.05,]

GO_overlap <- goana(as.vector(unique(overlap$target_entrez)))
KEGG_overlap <- kegga(as.vector(unique(overlap$target_entrez)))
GO_p_overlap <- GO_overlap[GO_overlap$P.DE < 0.05,]
KEGG_p_overlap <- KEGG_overlap[KEGG_overlap$P.DE < 0.05,]


#Determine miRNA genes and associated with osteoblast differentiation and ossification
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO term
attr<-listAttributes(ensembl)
filt <- listFilters(ensembl)
ossification.genes <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', "entrezgene_id"),
                            filters = 'go', values = 'GO:0043932', mart = ensembl)
obdiff.genes <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', "entrezgene_id"),
                      filters = 'go', values = 'GO:0001649', mart = ensembl)
reg_obdiff.genes <-getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', "entrezgene_id"),
                         filters = 'go', values = 'GO:0045667', mart = ensembl)
oss.genes <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', "entrezgene_id"),
                   filters = 'go', values = 'GO:0001503', mart = ensembl)

# Filter for targets 
ADRB2_3<- ADRB2_genesdfuni %>% filter(., target_entrez %in% ADRB2_table3$Var1)
goi <- as.character(c(ossification.genes$entrezgene_id,obdiff.genes$entrezgene_id, reg_obdiff.genes, oss.genes$entrezgene_id)) 
targets3<- ADRB2_3 %>% filter(., target_entrez %in% goi)
all_targets<- ADRB2_genesdfuni %>% filter(., target_entrez %in% goi)
