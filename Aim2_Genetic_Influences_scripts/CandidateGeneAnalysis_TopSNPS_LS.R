# Sensitivity Analysis LS
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
    lmmkin1 <- lmekin(formula = as.formula(paste("s8cbl24bd~gt_DS*BB +", cov, "+ (1|colnames(kmat))")), data = data$data[[i]],
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

# Extract Top SNPs ####
data_fem_nest_top_SNPs <- data_fem %>% filter(., POS %in% c(115803375, 240223080)) %>% group_by(., POS, Gene) %>% nest()
data_male_nest_top_SNPs <- data_male %>% filter(., POS %in% c(60001153, 60026732)) %>% group_by(., POS, Gene) %>% nest()

# Run Sex stratified model 
topSNP_results_f <- results_lmekin(data = data_fem_nest_top_SNPs, kmat = female_kmat, cov = "AGE8 + HGT8 + BMI8 + EST8")
int_topSNP_results_f <- topSNP_results_f %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)
topSNP_results_m <- results_lmekin(data = data_male_nest_top_SNPs, kmat = male_kmat, cov = "AGE8 + HGT8 + BMI8" )
int_topSNP_results_m <- topSNP_results_m %>% filter(., var == "gt_DS:BBYes") %>% filter(., p < 0.05)

cor.test(data1$s8cbl24bd, data1$f8cbnbmd)
table(is.na(data1$s8cbl24bd))

write.csv(topSNP_results_m, "TopSNPs_male_LS.csv", row.names = F, quote = F)
write.csv(topSNP_results_f, "TopSNPs_female_LS.csv", row.names = F, quote = F)
