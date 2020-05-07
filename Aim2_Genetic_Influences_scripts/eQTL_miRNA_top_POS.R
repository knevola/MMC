rm(list = ls())
library(tidyverse)
library(haven)
library(kinship2)
library(coxme)

# Read in Data ####
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/ktn_target_genes_4_2020")
topsnps<-read.csv("Significant_SNPs_GCTA_COJO.csv", stringsAsFactors = F)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
pedigree1 <- read.csv("share_ped_052517.csv")
load("CandidateGenes_4_2020.RData")

setwd("/home/clary@mmcf.mehealth.org/Framingham/eQTL Analysis/data")
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)
Colsums<- read.csv("ColumnSums.csv", stringsAsFactors = F)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

# Functions ####
extract_lmekin <- function(x, ...) {
  fcoef <- x$coefficients$fixed
  if (length(fcoef) > 0)  { # Not a ~1 model
    se <- sqrt(diag(x$var))
    tmp <- cbind(fcoef, se, round(fcoef/se,2),
                 signif(1 - pchisq((fcoef/ se)^2, 1), 2))
    dimnames(tmp) <- list(names(fcoef), c("Value", "Std Error", "z", "p")) 
  }
  results = data.frame(n = x$n, tmp)
}

results_miRNA_lmekin <- function(miRNA,data, kmat, cov){
  # miRNA - vector of miRNA names
  # data - all data, miRNA, top POS
  # kmat - kinship matrix
  # cov - vector of covariates
  results <- NULL
  for (j in 1:length(miRNA)){
    print(miRNA[j])
    results1 <- NULL
    for (i in 1:length(data$POS)){
      print(i)
      x <- Sys.time()
      lmmkin1 <- lmekin(formula = as.formula(paste(miRNA[j],"~gt_DS*BB+", cov, "+ (1|colnames(kmat))")), data = data$data[[i]],
                        varlist = list(kmat))
      new_results <- data.frame(extract_lmekin(lmmkin1), POS = data$POS[[i]], Gene = data$Gene[[i]], miRNA = miRNA[j])
      new_results$var <- row.names(new_results)
      results1 <- rbind(results1, new_results)
      y<-Sys.time()
      print(y-x)
    }
    results <- rbind(results, results1)
  }
  return(results)
}

# Merge in miRNA data ####
data<-pheno_snp_pos %>% unnest(., cols = c(data))
topdata <- data[data$POS %in% topsnps$bp,]


# Format miRNA technical variables
miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)
quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)

# Merge miRNA and top SNPs data
miRNA_pheno_tech <- merge(topdata, mirnatech, by.x = "Indiv", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno_tech, miRNAdat, by.x = "Indiv", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)
drop <- c('cvdpair', 'casecontrol', 'idtype.x', "idtype.y")
miRNA_pheno <- miRNA_pheno[, !(names(miRNA_pheno) %in% drop)]

# Set up Pedigree file for all subjects ####
data <- miRNA_pheno
missing_id <- as.numeric(setdiff(data$Indiv, pedigree1$shareid))
data1 <- pheno_snp_pos$data[[1]]
missing_sex <- data1$SEX[data1$Indiv %in% missing_id]
missing_sex[missing_sex == "Female"] <- 2
missing_sex[missing_sex == "Male"] <- 1
missing_sex <- as.numeric(missing_sex)
missing_ped <- data.frame(pedno = 1600:(1600+length(missing_id)-1), shareid = missing_id, fshare = NA, mshare = NA, sex = missing_sex, itwin = NA)
ped_all <- rbind(pedigree1, missing_ped)
ped<-with(ped_all, pedigree(id = shareid, fshare, mshare, sex=sex, famid=pedno))
kmat <- kinship(ped)

# Create sex stratified data and kmat ####
data_fem <- data %>% filter(., SEX == "Female")
data_male <- data %>% filter(., SEX == "Male")

kmat1<-as.matrix(kmat)
ids <- colnames(kmat1) %in% data1$Indiv
pheno_kmat <- kmat1[ids,ids]

female_ids <- colnames(kmat1) %in% data_fem$Indiv
female_kmat <- kmat1[female_ids,female_ids]

male_ids <- colnames(kmat1) %in% data_male$Indiv
male_kmat <- kmat1[male_ids,male_ids]

data_fem_nest <- data_fem %>% group_by(., POS, Gene) %>% nest()
data_male_nest <- data_male %>% group_by(., POS, Gene) %>% nest()

# Split miRNA into model types ####
Colsums1 <- Colsums[Colsums$NAs < 0.1 & !is.na(Colsums$NAs),] 
Colsums1 <- Colsums1 %>% filter(.,variables != 'idtype.y') %>% filter(.,variables != 'casecontrol')
Colsums2 <- Colsums[(Colsums$NAs > 0.9) & (Colsums$NAs < 0.95) & !is.na(Colsums$NAs),]
Colsums3 <- Colsums[(Colsums$NAs > 0.1 & Colsums$NAs < 0.9)& !is.na(Colsums$NAs),]
Colsums3 <- Colsums3 %>% filter(.,variables != 'cvdpair')

# Run top miRNA ~ SNPs, Sex-stratified models ####
fem_cov <- "AGE8 + HGT8 + BMI8 + EST8 + rankcon + rankqual + rank260"
male_cov <- "AGE8 + HGT8 + BMI8 + rankcon + rankqual + rank260"
top_miRNAs <- c("miR_19a_3p", "miR_186_5p_a1", "miR_186_5p_a2")
fem_results <- results_miRNA_lmekin(miRNA = top_miRNAs, data = data_fem_nest,kmat = female_kmat, cov = fem_cov)
male_results <- results_miRNA_lmekin(miRNA = top_miRNAs, data = data_male_nest,kmat = male_kmat, cov = male_cov)

# Filter results
fem_top_miR_sig <- fem_results %>% filter(., var %in% c("gt_DS", "gt_DS:BBYes")) %>% filter(.,p < 0.05 )
male_top_miR_sig <- male_results %>% filter(., var %in% c("gt_DS", "gt_DS:BBYes")) %>% filter(.,p < 0.05 )

# Run model 1 miRNA ~ SNPs, Sex Stratified Models ####
fem_results <- results_miRNA_lmekin(miRNA = Colsums1$X, data = data_fem_nest,kmat = female_kmat, cov = fem_cov)
male_results <- results_miRNA_lmekin(miRNA = Colsums1$X, data = data_male_nest,kmat = male_kmat, cov = male_cov)
fem_miRs_sig <- fem_results %>% filter(., var =="gt_DS:BBYes") %>% filter(.,p < 0.05 )

write.csv(fem_results,"Model1_miRNA_LMEKIN_eQTL_topSNPs_fem.csv", row.names = F, quote = F)
write.csv(male_results, "Model1_miRNA_LMEKIN_eQTL_topSNPs_male.csv", row.names = F, quote = F)


library(ggplot2)
ggplot(data_male_nest$data[[1]], aes(x = gt_DS, y = miR_186_5p_a2, col = BB)) + geom_line()
