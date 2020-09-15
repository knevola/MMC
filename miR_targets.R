rm(list = ls())
library(RmiR.Hs.miRNA)
library(RSQLite)
library(org.Hs.eg.db)

RSQLite::dbListTables(RmiR.Hs.miRNA_dbconn())
miRNAs <- c("hsa-miR-19a", "hsa-miR-184", "hsa-miR-222", "hsa-miR-125b")
# Tarbase: verified targets
tarbase <- dbReadTable(RmiR.Hs.miRNA_dbconn(), 'tarbase')[,1:2]
# TargetScan: next smallest database
targetscan <- dbReadTable(RmiR.Hs.miRNA_dbconn(), "targetscan")[,1:2]

miRNA_tarbase <- unique(tarbase[tarbase$mature_miRNA %in% miRNAs,]) # 4 genes
miRNA_targetscan <- unique(targetscan[targetscan$mature_miRNA %in% miRNAs,]) #1869 genes

hs <- org.Hs.eg.db
miR_TB_genes <- select(hs,keys = miRNA_tarbase$gene_id, 
                           columns = c("ENTREZID", "SYMBOL", 'GENENAME'),
                           keytype = "ENTREZID")
miR_TS_genes <- select(hs,keys = as.character(miRNA_targetscan$gene_id), 
                       columns = c("ENTREZID", "SYMBOL", 'GENENAME'),
                       keytype = "ENTREZID")
