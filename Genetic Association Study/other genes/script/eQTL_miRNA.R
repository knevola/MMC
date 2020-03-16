rm(list = ls())
library(tidyverse)
library(haven)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

topsnps<-read.csv("TopSNPsfinallist.csv", stringsAsFactors = F)
load("allSNPs_cleaned_andQC_genomewide.RData")

toppos <- topsnps$SNP
toppos <- as.integer(gsub("^..?:", "", toppos))
alldata <- tidy_clean2 %>% unnest(cols = c(data))
topdata <- alldata %>% filter(., POS %in% toppos) 

setwd("/home/clary@mmcf.mehealth.org/Framingham/eQTL Analysis/data")
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)
Colsums<- read.csv("ColumnSums.csv", stringsAsFactors = F)

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


# MLR Function miRNA~SNP
mlr_snp <- function(positions, miRNAs, snp_var, cov, data){
  mlr1 <- data.frame(NULL)
  microRNAs_col <- NULL
  position_col <- NULL
  for (i in 1:length(positions)){
    print("NEW POSITION")
    print(i)
    for(j in 1:length(miRNAs)){
      x <- lm(as.formula(paste(miRNAs[j],snp_var,cov, sep = "")), data = data[[2]][[i]])
      y <- summary(x)
      z<- as.data.frame(y$coefficients)
      z$variable <- row.names(z)
      z$FDR <- p.adjust(z$`Pr(>|t|)`, method = "BH", n =  158)
      row.names(z)<- NULL
      mlr1 <- rbind(mlr1, z)
      microRNAs_col <- c(microRNAs_col,rep(miRNAs[j], nrow(z)))
      position_col <- c(position_col,rep(positions[i], nrow(z)))
      
    }
  }
  mlr1$miRNA <- microRNAs_col
  mlr1$position <- position_col
  return(mlr1)
}

# MLR Function (Binary) for Model 2 and 3
mlr_snp2_3 <- function(positions, miRNAs, snp_var, cov, data){
  mlr1 <- data.frame(NULL)
  microRNAs_col <- NULL
  position_col <- NULL
  for (i in 1:length(positions)){
    print("NEW POSITION")
    print(i)
    for(j in 1:length(miRNAs)){
      x <- glm(as.formula(paste(miRNAs[j],snp_var,cov, sep = "")), data = data[[2]][[i]])
      y <- summary(x)
      z<- as.data.frame(y$coefficients)
      z$variable <- row.names(z)
      row.names(z)<- NULL
      mlr1 <- rbind(mlr1, z)
      microRNAs_col <- c(microRNAs_col,rep(miRNAs[j], nrow(z)))
      position_col <- c(position_col,rep(positions[i], nrow(z)))
    }
  }
  mlr1$miRNA <- microRNAs_col
  mlr1$position <- position_col
  return(mlr1)
}

#Fisher p-value calculation function for model 3 (Need to fix!)####
snp_fisher <- function(data_lin, data_log){
  # Combine pos and miRNA columns
  data_lin$Combined <- paste(data_lin$miRNA, data_lin$position, data_lin$variable, sep = " : ")
  data_log$Combined <- paste(gsub(x = data_log$miRNA,pattern = "_nas",replacement = ''), data_log$position,
                             data_log$variable,sep = " : ")
  
  # Calculate Fisher p-values####
  names(data_lin) <- c("Estimate", "SE", "tvalue", "pvalue", "variable","miRNA", "position", "FDR", "Combined")
  names(data_log) <- c("Estimate", "SE", "tvalue", "pvalue", "variable","miRNA", "position", "FDR", "Combined")
  model1a2<- merge(data_lin, data_log, by = "Combined")
  keep1 <- c("pvalue.x","pvalue.y") 
  model1a26<- model1a2[names(model1a2) %in% keep1]
  results<-fisher.method(model1a26, p.corr = "BH")
  model3<-merge(model1a2,results, by.x = 0, by.y = 0)
  return(model3)
}

# Split miRNA into model types 
Colsums1 <- Colsums[Colsums$NAs < 0.1 & !is.na(Colsums$NAs),] 
Colsums1 <- Colsums1 %>% filter(.,variables != 'idtype.y') %>% filter(.,variables != 'casecontrol')
Colsums2 <- Colsums[(Colsums$NAs > 0.9) & (Colsums$NAs < 0.95) & !is.na(Colsums$NAs),]
Colsums3 <- Colsums[(Colsums$NAs > 0.1 & Colsums$NAs < 0.9)& !is.na(Colsums$NAs),]
Colsums3 <- Colsums3 %>% filter(.,variables != 'cvdpair')

# Group data by position
nested_data <- miRNA_pheno %>% group_by(., POS) %>% nest()

# Covariate names/input for function
positions <- as.character(nested_data[[1]])
pheno_cov <- '+ SEX + AGE8 + HGT8 + WGT8'
tech_cov <- ' + rankcon + rankqual + rank260'

# Run Model 1 Best Guess ####
mlr1<-mlr_snp(positions = positions, miRNAs = Colsums1$variables, snp_var = '~ Genotype', cov = paste(pheno_cov, tech_cov), data = nested_data)
miRs <- mlr1 %>%  
  filter(., variable %in% c("Genotypeheterozygous", "Genotypehomozygous alternative")) %>% filter(., FDR < 0.05)
mirs5 <- mlr1 %>%  
  filter(., variable %in% c("Genotypeheterozygous", "Genotypehomozygous alternative")) %>% 
  filter(., `Pr(>|t|)` < 0.05)

# Position mapping
pos <- topsnps %>% select(., Gene, SNP)
pos$POS <- gsub(pattern = "^..?:", replacement = "", x = pos$SNP)
miRs <- merge(miRs, pos, by.x = "position", by.y = "POS")
mirs5 <- merge(mirs5, pos, by.x = "position", by.y = "POS")
write.csv(miRs, "FDRsig_eqtl.csv", row.names = F, quote = F)
write.csv(mirs5,"psig_eqtl.csv", row.names = F, quote = F)

df<- as.data.frame(table(miRs$Gene,miRs$miRNA))
df1 <- df[df$Freq != 0,]
write.csv(df1, "SignificantMiRNAswithGenes.csv")
