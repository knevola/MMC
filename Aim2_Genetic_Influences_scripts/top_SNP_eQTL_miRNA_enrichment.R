# Top SNPs eQTL miRNA enrichment analysis
rm(list = ls())
library(tidyverse)
library(multiMiR)
library(limma)
library(ggplot2)
library(biomaRt)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

# Read files
sig_mirs<-read.csv("Female_Top2_SNPs_eQTL_sig.csv")

# Format miRNAs to be read by multimiR
sig_mirs$miRNA <- gsub(pattern = "miR", replacement = "hsa-miR",x = sig_mirs$miRNA)
sig_mirs$miRNA <- gsub(pattern = "_", replacement = "-", x=sig_mirs$miRNA)

# Separate into SNPs
ADRB1_mirs <- sig_mirs %>% filter(., sig_mirs$Gene == "ADRB1")
HDAC4_mirs <- sig_mirs %>% filter(., sig_mirs$Gene == "HDAC4")

# Determine targets
ADRB1_miR_targets<- get_multimir(org = "hsa",mirna = ADRB1_mirs$miRNA)
ADRB1_miR_targets <- ADRB1_miR_targets@data
ADRB1_miR_targets <- unique(ADRB1_miR_targets[3:5])

HDAC4_miR_targets<- get_multimir(org = "hsa",mirna = HDAC4_mirs$miRNA)
HDAC4_miR_targets <- HDAC4_miR_targets@data
HDAC4_miR_targets <- unique(HDAC4_miR_targets[3:5])

# Filter for targetted by more that 1 miRNA
ADRB1_target_freq<-as.data.frame(table(ADRB1_miR_targets$target_entrez))
ADRB1_target_freq_3 <- ADRB1_target_freq[(ADRB1_target_freq$Freq > 2) & (ADRB1_target_freq$Freq < 7),]
ADRB1_target_freq_3$Var1 <- as.character(ADRB1_target_freq_3$Var1)

HDAC4_target_freq<-as.data.frame(table(HDAC4_miR_targets$target_entrez))
HDAC4_target_freq_3 <- HDAC4_target_freq[(HDAC4_target_freq$Freq > 2) & (HDAC4_target_freq$Freq < 7),]
HDAC4_target_freq_3$Var1 <- as.character(HDAC4_target_freq_3$Var1)

# Enrichment analysis
ADRB1_GO <- goana(de = unique(ADRB1_target_freq_3$Var1))
ADRB1_GO_sig <- ADRB1_GO %>% filter(., P.DE < 0.05)
ADRB1_KEGG <- kegga(de = ADRB1_target_freq_3$Var1)
ADRB1_KEGG_sig <- ADRB1_KEGG %>% filter(., P.DE < 0.05)

HDAC4_GO <- goana(de = unique(HDAC4_target_freq_3$Var1))
HDAC4_GO_sig <- HDAC4_GO %>% filter(., P.DE < 0.05)
HDAC4_KEGG <- kegga(de = HDAC4_target_freq_3$Var1)
HDAC4_KEGG_sig <- HDAC4_KEGG %>% filter(., P.DE < 0.05)

# Filter osteo or bone related terms
bone_ADRB1_GO<- ADRB1_GO_sig %>% filter(., Ont == "BP")%>%  filter(.,str_detect(Term, "osteo*")& !str_detect(Term, "cortico")|str_detect(Term,"bone")|str_detect(Term,"estrogen")|str_detect(Term,"insulin")|str_detect(Term,"thyroid"))
bone_ADRB1_KEGG<- ADRB1_KEGG_sig %>%  filter(.,str_detect(Pathway, "Osteo*")& !str_detect(Pathway, "Cortico")|str_detect(Pathway,"Bone")|str_detect(Pathway,"Estrogen")|str_detect(Pathway,"Insulin")|str_detect(Pathway,"Thyroid"))
bone_HDAC4_GO<- HDAC4_GO_sig %>% filter(., Ont == "BP")%>%  filter(.,str_detect(Term, "osteo*")& !str_detect(Term, "cortico")|str_detect(Term,"bone")|str_detect(Term,"estrogen")|str_detect(Term,"insulin")|str_detect(Term,"thyroid"))
bone_HDAC4_KEGG<- HDAC4_KEGG_sig %>%  filter(.,str_detect(Pathway, "Osteo*")& !str_detect(Pathway, "Cortico")|str_detect(Pathway,"Bone")|str_detect(Pathway,"Estrogen")|str_detect(Pathway,"Insulin")|str_detect(Pathway,"Thyroid"))

bone_ADRB1_GO$Gene <- "ADRB1"
bone_ADRB1_GO$SNP <- "rs12414657"
bone_HDAC4_GO$Gene <- "HDAC4"
bone_HDAC4_GO$SNP <- "rs11124190"
bone_GO <- rbind(bone_ADRB1_GO, bone_HDAC4_GO)
write.csv(bone_GO, "miRNA_GO.csv", row.names = F, quote = F)

# Genes associated with annotations of interest
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, entrez_id and go_id for all genes annotated with GO term
z<-listFilters(ensembl)
insulin_sig <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'go_id'),
                   filters = 'go_parent_name', values = 'insulin receptor signaling pathway', mart = ensembl)
insulin_sig_genes <-insulin_sig[insulin_sig$entrezgene_id %in% ADRB1_target_freq_3$Var1,]
unique(insulin_sig_genes$hgnc_symbol)

Oc_sig <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'go_id'),
                     filters = 'go_parent_name', values = 'osteoclast differentiation', mart = ensembl)
Oc_sig_genes <-Oc_sig[Oc_sig$entrezgene_id %in% ADRB1_target_freq_3$Var1,]
unique(Oc_sig_genes$hgnc_symbol)

write.cs
