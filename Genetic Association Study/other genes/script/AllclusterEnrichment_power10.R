# Enrichment of all clusters mRNA and miRNA
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

library(dplyr)
library(multiMiR)
library(edgeR)
library(clusterProfiler)
library(readr)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
miRNA <- read.csv("MiRNA_significance_module_membership_nofiltergenetics10.csv")

miRtarget <- function(miRNAs){
  miRNAs <- gsub(x = miRNAs,pattern ="_",replacement = "-" )
  miRNAs <- gsub(x = miRNAs, pattern = "miR", replacement = "hsa-miR")
  miRNAs <- gsub(x = miRNAs, pattern = "5p-a2", replacement = "5p")
  miRNAs <- gsub(x = miRNAs, pattern = "5p-a1", replacement = "5p")
  genes<- get_multimir(mirna = miRNAs, summary = TRUE)
  genes <- genes@data
  return(genes)
}

miRfilter <- function(genes){
  table <- as.data.frame(table(genes$target_entrez))
  table2 <- table[table$Freq > 2,]
  genes2 <- genes %>% filter(., target_entrez %in% table2$Var1)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

GOmRNA <- function(mRNAlist, universe){
  GO <- goana(de = mRNAlist, universe = universe, FDR = 0.05)
  GO_p <- GO[GO$P.DE < 0.001,]
}

KEGGmRNA <- function(mRNAlist, universe){
  KEGG <- kegga(de = mRNAlist, universe = universe, FDR = 0.05)
  KEGG_p <- KEGG[KEGG$P.DE < 0.001,]
}

# MIR ENRICHMENT ####

# Separated in miRNA clusters
table(miRNA$mergedColors)
miRNA_blue <- miRNA %>% filter(., mergedColors == "blue")
miRNA_brown <- miRNA %>% filter(., mergedColors == "brown")
miRNA_grey <- miRNA %>% filter(., mergedColors == "grey")
miRNA_turquoise <- miRNA %>% filter(., mergedColors == "turquoise")
miRNA_yellow <- miRNA %>% filter(., mergedColors == "yellow")

# Gene Targets of Each Cluster
genetargets_blue <- miRtarget(miRNAs = as.character(miRNA_blue$miRNA))
genetargets_brown <- miRtarget(miRNAs = as.character(miRNA_brown$miRNA))
genetargets_grey <- miRtarget(miRNAs = as.character(miRNA_grey$miRNA))
genetargets_turquoise <- miRtarget(miRNAs = as.character(miRNA_turquoise$miRNA))
genetargets_yellow <- miRtarget(miRNAs = as.character(miRNA_yellow$miRNA))

# Gene Filter
genetargets_blue2 <- miRfilter(genetargets_blue)
genetargets_brown2 <- miRfilter(genetargets_brown)
genetargets_grey2 <- miRfilter(genetargets_grey)
genetargets_turquoise2 <- miRfilter(genetargets_turquoise)
genetargets_yellow2 <- miRfilter(genetargets_yellow)


#GO miRs
GO_blue <- goana(genetargets_blue2$target_entrez, FDR = 0.05)
GO_brown <- goana(genetargets_brown2$target_entrez, FDR = 0.05)
GO_grey <- goana(genetargets_grey2$target_entrez, FDR = 0.05)
GO_turquoise <- goana(genetargets_turquoise2$target_entrez, FDR = 0.05)
GO_yellow <- goana(genetargets_yellow2$target_entrez, FDR = 0.05)

# Filter GO
GO_blue_p <- GO_blue[GO_blue$P.DE < 0.05,]
GO_brown_p <- GO_brown[GO_brown$P.DE < 0.05,]
GO_grey_p <- GO_grey[GO_grey$P.DE < 0.05,]
GO_turquoise_p <- GO_turquoise[GO_turquoise$P.DE < 0.05,]
GO_yellow_p <- GO_yellow[GO_yellow$P.DE < 0.05,]

GO_blue_only<-GO_blue_p[GO_blue_p$Term %!in% c(GO_brown_p$Term, GO_grey_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_brown_only <- GO_brown_p[GO_brown_p$Term %!in% c(GO_blue_p$Term, GO_grey_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_grey_only <- GO_grey_p[GO_grey_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_turquoise_p$Term, GO_yellow_p$Term),]
GO_turquoise_only <- GO_turquoise_p[GO_turquoise_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_grey_p$Term, GO_yellow_p$Term),]
GO_yellow_only <- GO_yellow_p[GO_yellow_p$Term %!in% c(GO_blue_p$Term, GO_brown_p$Term, GO_grey_p$Term, GO_turquoise_p$Term),]

GO_blue_only1 <- GO_blue_only %>% filter(., P.DE < 0.001)
GO_brown_only1 <- GO_brown_only %>% filter(., P.DE < 0.001)
GO_turquoise_only1 <- GO_turquoise_only %>% filter(., P.DE < 0.001)
GO_yellow_only1 <- GO_yellow_only %>% filter(., P.DE < 0.001)

# KEGG miRs 
KEGG_blue <- kegga(genetargets_blue2$target_entrez)
KEGG_brown <- kegga(genetargets_brown2$target_entrez)
KEGG_grey <- kegga(genetargets_grey2$target_entrez)
KEGG_turquoise <- kegga(genetargets_turquoise2$target_entrez)
KEGG_yellow <- kegga(genetargets_yellow2$target_entrez)

# Filter KEGG
KEGG_blue_p <- KEGG_blue[KEGG_blue$P.DE < 0.05,]
KEGG_brown_p <- KEGG_brown[KEGG_brown$P.DE < 0.05,]
KEGG_grey_p <- KEGG_grey[KEGG_grey$P.DE < 0.05,]
KEGG_turquoise_p <- KEGG_turquoise[KEGG_turquoise$P.DE < 0.05,]
KEGG_yellow_p <- KEGG_yellow[KEGG_yellow$P.DE < 0.05,]

KEGG_blue_only <- KEGG_blue_p[KEGG_blue_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_grey_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_brown_only <- KEGG_brown_p[KEGG_brown_p$Pathway %!in% c(KEGG_blue_p$Pathway, KEGG_grey_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_grey_only <- KEGG_grey_p[KEGG_grey_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_blue_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_turquoise_only <- KEGG_turquoise_p[KEGG_turquoise_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_grey_p$Pathway, KEGG_blue_p$Pathway, KEGG_yellow_p$Pathway),]
KEGG_yellow_only <- KEGG_yellow_p[KEGG_yellow_p$Pathway %!in% c(KEGG_brown_p$Pathway, KEGG_grey_p$Pathway, KEGG_turquoise_p$Pathway, KEGG_blue_p$Pathway),]



