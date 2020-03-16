# Read in all .Rdata files
rm(list = ls())
library(tidyverse)
library(xlsx)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes")

genotype_files<- list.files("./", pattern = "_filtered_genotype.RData")
pheno <- read.csv("PhenoData_5_28.csv")

load(genotype_files[1])
adra1d_tidy <- tidy
adra1d  <- adra1d_tidy$gt
adra1d$Gene <- "ADRA1D"
length(unique(adra1d$POS))

load(genotype_files[2])
adrb1 <- tidy$gt
adrb1$Gene <- "ADRB1"
length(unique(adrb1$POS))

load(genotype_files[3])
adrb2_tidy <- tidy
adrb2 <- adrb2_tidy$gt
adrb2$Gene <- "ADRB2"
length(unique(adrb2$POS))


load(genotype_files[4])
bmp2_tidy <- tidy
bmp2 <- bmp2_tidy$gt
bmp2$Gene <- "BMP2"
length(unique(bmp2$POS))


load(genotype_files[5])
bmp4_tidy <- tidy
bmp4 <- bmp4_tidy$gt
bmp4$Gene <- "BMP4"
length(unique(bmp4$POS))

load(genotype_files[6])
bmp5_tidy <- tidy
bmp5 <- bmp5_tidy$gt
bmp5$Gene <- "BMP5"
length(unique(bmp5$POS))


load(genotype_files[7])
bmp6_tidy <- tidy
bmp6 <- bmp6_tidy$gt
bmp6$Gene <- "BMP6"
length(unique(bmp6$POS))

load(genotype_files[8])
creb5_tidy <- tidy
creb5 <- creb5_tidy$gt
creb5$Gene <- "CREB5"
length(unique(creb5$POS))

load(genotype_files[9])
crebbp_tidy <- tidy
crebbp <- crebbp_tidy$gt
crebbp$Gene <- "CREBBP"
length(unique(crebbp$POS))

load(genotype_files[10])
foxO1_tidy <- tidy
foxO1 <- foxO1_tidy$gt
foxO1$Gene <- "FoxO1"
length(unique(foxO1$POS))


load(genotype_files[11])
hdac4_tidy <- tidy
hdac4 <- hdac4_tidy$gt
hdac4$Gene <- "HDAC4"
length(unique(hdac4$POS))


load(genotype_files[12])
prkar2a_tidy <-tidy
prkar2a <- prkar2a_tidy$gt
prkar2a$Gene <- "PRKAR2A"
length(unique(prkar2a$POS))

load(genotype_files[13])
rank_tidy <- tidy
rank <- rank_tidy$gt
rank$Gene <- "RANK"
length(unique(rank$POS))


load(genotype_files[14])
opg_tidy <- tidy
opg <- opg_tidy$gt
opg$Gene <- "OPG"
length(unique(opg$POS))

load(genotype_files[15])
rankl_tidy <- tidy
rankl <- rankl_tidy$gt
rankl$Gene <- "RANKL"
length(unique(rankl$POS))


# Merge genotype files
all_dfs<- list(adra1d, adrb1, adrb2, bmp2, bmp4, bmp5, bmp6, creb5, crebbp, foxO1, hdac4, opg, prkar2a, rank, rankl)
all_genes <- bind_rows(all_dfs)

# Filter for people with phenotype data and merge with pheno data
all_genes_pheno <- merge(all_genes, pheno, by.x = "Indiv", by.y = "shareid")

# Number of SNPs (4078)
length(unique(all_genes_pheno$POS)) 

# Number of People (1527)
length(unique(all_genes_pheno$Indiv)) 

# Check that bad snps were filtered out
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
good_snps<-read.table("WellImputedPositions8.txt", header = T, sep = "\t")
bad_snps <- setdiff(unique(all_genes_pheno$POS),as.integer(good_snps$POS)) # Empty so no good snps made it this far

pheno_snp <- within(all_genes_pheno, Genotype <- relevel(Genotype, ref = "normal"))
pheno_snp_pos<- pheno_snp %>% group_by(.,POS) %>% nest()
save.Object(pheno_snp_pos, "All_Genes_SNPs.RData")
