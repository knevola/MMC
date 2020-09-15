rm(list= ls())
library(kinship2)
library(tidyverse)
library(coxme)
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
results_lmekin <- function(data, kmat, cov){
  results <- NULL
  for (i in 1:length(data$POS)){
    print(i)
    x <- Sys.time()
    lmmkin1 <- lmekin(formula = as.formula(paste("f8cbnbmd~gt_DS*BB +", cov, "+ (1|colnames(kmat))")), data = data$data[[i]],
                      varlist = list(kmat))
    new_results <- data.frame(extract_lmekin(lmmkin1), POS = data$POS[[i]], Gene = data$Gene[[i]])
    new_results$var <- row.names(new_results)
    results <- rbind(results, new_results)
    y<-Sys.time()
    print(y-x)
  }
  return(results)
}

# Read in Files ####
pedigree1 <- read.csv("share_ped_052517.csv")
load("CandidateGenes_4_2020.RData")

# Set up Pedigree file for all subjects
data<-pheno_snp_pos %>% unnest(., cols = c(data))

missing_id <- as.numeric(setdiff(data$Indiv, pedigree1$shareid))
data1 <- pheno_snp_pos$data[[1]]
missing_sex <- data1$SEX[data1$Indiv %in% missing_id]
missing_sex[missing_sex == "Female"] <- 2
missing_sex[missing_sex == "Male"] <- 1
missing_sex <- as.numeric(missing_sex)
missing_ped <- data.frame(pedno = 1600:1680, shareid = missing_id, fshare = NA, mshare = NA, sex = missing_sex, itwin = NA)
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

#  Are there any Identical twins? ####
pheno_pedigree <- ped_all[ped_all$shareid %in% data1$Indiv,]
twins <- pheno_pedigree[is.na(pheno_pedigree$itwin)== F,] # No sets of identical twin in the dataset

# ADRB1 and ADRB2 ####
# Set Data
data_fem_nest_ADRB1_ADRB2 <- data_fem %>% filter(., Gene %in% c("ADRB1", "ADRB2")) %>% group_by(., POS, Gene) %>% nest()
data_male_nest_ADRB1_ADRB2 <- data_male %>% filter(., Gene %in% c("ADRB1", "ADRB2")) %>% group_by(., POS, Gene) %>% nest()
data_nest_ADRB1_ADRB2 <- data %>% filter(., Gene %in% c("ADRB1", "ADRB2")) %>% group_by(., POS, Gene) %>% nest()

# Run Sex stratified model 
adrb1_adrb2_results_f <- results_lmekin(data = data_fem_nest_ADRB1_ADRB2, kmat = female_kmat, cov = "AGE8 + HGT8 + BMI8 + EST8")
int_adrb1_adrb2_results_f <- adrb1_adrb2_results_f %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)
adrb1_adrb2_results_m <- results_lmekin(data = data_male_nest_ADRB1_ADRB2, kmat = male_kmat, cov = "AGE8 + HGT8 + BMI8" )
int_adrb1_adrb2_results_m <- adrb1_adrb2_results_m %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)

# Run sex-combined model
adrb1_adrb2_results <- results_lmekin(data = data_nest_ADRB1_ADRB2, kmat = pheno_kmat, cov = "AGE8 + HGT8 + BMI8 + SEX" )
int_adrb1_adrb2_results <- adrb1_adrb2_results %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)

# All SNPs Female Model ####
data_fem_nest <- data_fem %>% group_by(., POS, Gene) %>% nest()
fem_results <- results_lmekin(data = data_fem_nest, kmat = female_kmat, cov = "AGE8 + HGT8 + BMI8 + EST8")
int_fem_results <- fem_results %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)
  
write.csv(fem_results, "AllSNPs_LMEKIN_results_female.csv", quote = F, row.names = F)
write.csv(int_fem_results, "AllSNP_sig_interaction_LMEKIN_female.csv", quote = F, row.names = F)
table(int_fem_results$Gene)

int_fem_results_15 <- fem_results %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.003333333)
table(int_fem_results_15$Gene)

#All SNPs Male Model
data_male_nest <- data_male %>% group_by(., POS, Gene) %>% nest()
male_results <- results_lmekin(data = data_male_nest, kmat = male_kmat, cov = "AGE8 + HGT8 + BMI8")
int_male_results <- male_results %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)

int_male_results_15 <- male_results %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.003333333)

write.csv(male_results, "AllSNPs_LMEKIN_results_male.csv", quote = F, row.names = F)
write.csv(int_male_results, "AllSNP_sig_interaction_LMEKIN_male.csv", quote = F, row.names = F)

