# SNP Analysis
rm(list = ls())
library(haven)
library(dplyr)
library(metaseqR)
library(VennDiagram)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv("/home/clary@mmcf.mehealth.org/Framingham/HEglinton/PhenoData_wsnp.csv")
# X5.148206646
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

# Merge pheno with miRNA ####
quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)
Colsums <- data.frame(variables = names(miRNA_pheno), NAs = colSums(is.na(miRNA_pheno))/nrow(miRNA_pheno))
Colsums$NAs[1:34]<- NA 
drop <- c('cvdpair', 'casecontrol', 'idtype.x', "idtype.y")
miRNA_pheno <- miRNA_pheno[, !(names(miRNA_pheno) %in% drop)]
miRNA_pheno$X5.148206646 <- as.factor(miRNA_pheno$X5.148206646)
miRNA_pheno$X10.115805056 <- as.factor(miRNA_pheno$X10.115805056)
# Setup model specific data ####
# Model 1 Only (linear)
Colsums1 <- Colsums[Colsums$NAs < 0.1 | is.na(Colsums$NAs) == T,]
miRNA_pheno1 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums1$variables] # Filter out miRNA with more than 10% NAs
names1 <- names(miRNA_pheno1)[-c(1:33)]
# Model 2 only (log/binary)
Colsums2 <- Colsums[Colsums$NAs > 0.9 | is.na(Colsums$NAs) == T,]
miRNA_pheno2 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums2$variables]
miRNAdatna <- miRNA_pheno2[-c(1:33)]
miRNAdatna[!is.na(miRNAdatna)] <- 0 #Detectable
miRNAdatna[is.na(miRNAdatna)] <- 1 #Not Detectable
miRNAdatna[sapply(miRNAdatna, is.numeric)] <- lapply(miRNAdatna[sapply(miRNAdatna, is.numeric)], as.logical)
miRNA_pheno2 <- cbind(miRNA_pheno2[c(1:33)], miRNAdatna)
names2 <- names(miRNA_pheno2)[-c(1:33)]

# Model 3:Both Model 1 and 2
Colsums3 <- Colsums[(Colsums$NAs > 0.1 & Colsums$NAs < 0.9)| is.na(Colsums$NAs) == T,]
miRNA_pheno3 <- miRNA_pheno[names(miRNA_pheno) %in% Colsums3$variables]
miRNAdatna <- miRNA_pheno3[-c(1:33)]
miRNAdatna[!is.na(miRNAdatna)] <- 0 #Detectable
miRNAdatna[is.na(miRNAdatna)] <- 1 #Not Detectable
miRNAdatna[sapply(miRNAdatna, is.numeric)] <- lapply(miRNAdatna[sapply(miRNAdatna, is.numeric)], as.logical)
miRNA_pheno3_2 <- cbind(miRNA_pheno3[c(1:33)], miRNAdatna)
names3 <- names(miRNA_pheno3)[-c(1:33)]

# MLR Function ####
mlr<- function(names, cov, data, var){
  name <- NULL
  mlr1 <- NULL
  for (i in 1:length(names)){
    #tryCatch({
    print(i)
    x <- lm(as.formula(paste(names[i], cov, sep = "")), data = data)
    y <- summary(x)
    z<- as.data.frame(y$coefficients)
    z$variable <- row.names(z)
    row.names(z)<- NULL
    if (i == 1){
      mlr1 <- z
      name <- rep(names[i], nrow(z))
    }else{
      mlr1 <- rbind(mlr1, z)
      name <- c(name, rep(names[i], nrow(z)))
    }
    #}, error = function(e){})
  }
  mlr_results <- cbind(miRNA = name, mlr1)
  mlr_fdr <- mlr_results %>% filter(., variable == var)
  mlr_fdr$FDR <- p.adjust(p = mlr_fdr$`Pr(>|t|)`, method = "BH", n = length(names)) 
  mlr_results <- merge(x = mlr_results, y = mlr_fdr, all.x = T)
  return(mlr_results)
}

fisherpvalue <- function(data1, data2, var){
  names(data1) <- c("miRNA", "Estimate", "SE", "tvalue", "pvalue", "variable", "FDR")
  names(data2) <- c("miRNA", "Estimate", "SE", "tvalue", "pvalue", "variable", "FDR")
  data1 <- data1 %>% filter(.,variable == var)
  data2 <- data2 %>%  filter(.,variable == var)
  model1a2<- merge(data1, data2, by.x = 'miRNA', by.y = "miRNA")
  row.names(model1a2)<- model1a2$miRNA
  keep1 <- c("pvalue.x","pvalue.y") 
  model1a26<- model1a2[names(model1a2) %in% keep1]
  results<-fisher.method(model1a26, p.corr = "BH")
  model3<-merge(model1a2,results, by.x = 0, by.y = 0)
  return(model3)
}

SNPcov <-"~X5.148206646 + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"
SNPBBcov <- "~X5.148206646 * BB + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"
SNP2cov <- "~X10.115805056 + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"
SNP2BBcov <- "~X10.115805056*BB + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260"

SNPmlr1_1 <- mlr(names = names1, cov = SNPcov, data = miRNA_pheno1, var = "X5.1482066461")
SNPmlr1_2 <- mlr(names = names1, cov = SNPcov, data = miRNA_pheno1, var = "X5.1482066462")
SNPBBmlr1_1 <- mlr(names = names1, cov = SNPBBcov, data = miRNA_pheno1, var = "X5.1482066461:BBYes" )
SNPBBmlr1_2 <- mlr(names = names1, cov = SNPBBcov, data = miRNA_pheno1, var = "X5.1482066462:BBYes" )

SNPmlr1_1_p <- SNPmlr1_1 %>% filter(., variable %in% c("X5.1482066461", "X5.1482066462") )  %>% filter(., `Pr(>|t|)` < 0.05)
SNPBBmlr1_1_p <- SNPBBmlr1_1 %>% filter(., variable %in% c("X5.1482066461:BBYes", "X5.1482066462:BBYes") )  %>% filter(., `Pr(>|t|)` < 0.05)

SNPmlr2_1 <- mlr(names = names2, cov = SNPcov, data = miRNA_pheno2, var = "X5.1482066461")
SNPmlr2_2 <- mlr(names = names2, cov = SNPcov, data = miRNA_pheno2, var = "X5.1482066462")
SNPBBmlr2_1 <- mlr(names = names2, cov = SNPBBcov, data = miRNA_pheno2, var = "X5.1482066461:BBYes" )
SNPBBmlr2_2 <- mlr(names = names2, cov = SNPBBcov, data = miRNA_pheno2, var = "X5.1482066462:BBYes" )

SNPmlr2_1_p <- SNPmlr2_1 %>% filter(., variable %in% c("X5.1482066461", "X5.1482066462") )  %>% filter(., `Pr(>|t|)` < 0.05)
SNPBBmlr2_1_p <- SNPBBmlr2_1 %>% filter(., variable %in% c("X5.1482066461:BBYes", "X5.1482066462:BBYes") )  %>% filter(., `Pr(>|t|)` < 0.05)


SNPmlr3_1 <- mlr(names = names3, cov = SNPcov, data = miRNA_pheno3, var = "X5.1482066461")
SNPmlr3_2 <- mlr(names = names3, cov = SNPcov, data = miRNA_pheno3_2,var = "X5.1482066461" )
SNPmlr3_x1 <-fisherpvalue(data1 = SNPmlr3_1, data2 = SNPmlr3_2, var = "X5.1482066461")
SNPmlr3_x2 <-fisherpvalue(data1 = SNPmlr3_1, data2 = SNPmlr3_2, var = "X5.1482066462")

SNPBBmlr3_1 <- mlr(names = names3, cov = SNPBBcov, data = miRNA_pheno3, var = "X5.1482066461:BBYes" )
SNPBBmlr3_2 <- mlr(names = names3, cov = SNPBBcov, data = miRNA_pheno3_2, var = "X5.1482066461:BBYes" )
SNPBBmlr3_x1 <-fisherpvalue(data1 = SNPBBmlr3_1, data2 = SNPBBmlr3_2, var = "X5.1482066461:BBYes")
SNPBBmlr3_x2 <-fisherpvalue(data1 = SNPBBmlr3_1, data2 = SNPBBmlr3_2, var = "X5.1482066462:BBYes")

SNPmlr3_1_p <- SNPmlr3_x1 %>%  filter(.,p.value < 0.05)
SNPmlr3_2_p <- SNPmlr3_x2 %>%  filter(.,p.value < 0.05)
SNPBBmlr3_1_p <- SNPBBmlr3_x1 %>% filter(.,p.value < 0.05)
SNPBBmlr3_2_p <- SNPBBmlr3_x2 %>% filter(.,p.value < 0.05)

# No clinical variables ####
SNPcov1 <- "~X5.148206646"
SNPBBcov1 <- "~X5.148206646 * BB"

SNP_clean1_1 <- mlr(names = names1,cov = SNPcov1, data = miRNA_pheno1, var = "X5.1482066461")
SNP_clean1_2 <- mlr(names = names1,cov = SNPcov1, data = miRNA_pheno1, var = "X5.1482066462")
SNPBB_clean1_1 <- mlr(names = names1,cov = SNPBBcov1, data = miRNA_pheno1, var = "X5.1482066461:BBYes")
SNPBB_clean1_2 <- mlr(names = names1,cov = SNPBBcov1, data = miRNA_pheno1, var = "X5.1482066462:BBYes")

SNP_clean1_1_p <- SNP_clean1_1 %>% filter(., variable %in% c("X5.1482066461", "X5.1482066462") )  %>% filter(., `Pr(>|t|)` < 0.05)
SNPBB_clean1_1_p <- SNPBB_clean1_1 %>% filter(., variable %in% c("X5.1482066461:BBYes", "X5.1482066462:BBYes") )  %>% filter(., `Pr(>|t|)` < 0.05)
intersect(SNP_clean1_1_p$miRNA, SNPBB_clean1_1_p$miRNA)

SNP_clean2_1 <- mlr(names = names2,cov = SNPcov1, data = miRNA_pheno2, var = "X5.1482066461")
SNP_clean2_2 <- mlr(names = names2,cov = SNPcov1, data = miRNA_pheno2, var = "X5.1482066462")
SNPBB_clean2_1 <- mlr(names = names2,cov = SNPBBcov1, data = miRNA_pheno2, var = "X5.1482066461:BBYes")
SNPBB_clean2_2 <- mlr(names = names2,cov = SNPBBcov1, data = miRNA_pheno2, var = "X5.1482066462:BBYes")

SNP_clean2_2_f <- SNP_clean2_2 %>% filter(., FDR < 0.1)
SNPBB_clean2_2_f <- SNPBB_clean2_2 %>% filter(., FDR < 0.1)
SNP_clean2_1_p <- SNP_clean2_1 %>% filter(., variable %in% c("X5.1482066461", "X5.1482066462") )  %>% filter(., `Pr(>|t|)` < 0.05)
SNPBB_clean2_1_p <- SNPBB_clean2_1 %>% filter(., variable %in% c("X5.1482066461:BBYes", "X5.1482066462:BBYes") )  %>% filter(., `Pr(>|t|)` < 0.05)

SNPclean3_1 <- mlr(names = names3, cov = SNPcov1, data = miRNA_pheno3, var = "X5.1482066461")
SNPclean3_2 <- mlr(names = names3, cov = SNPcov1, data = miRNA_pheno3_2,var = "X5.1482066461" )
SNPclean3_x1 <-fisherpvalue(data1 = SNPclean3_1, data2 = SNPclean3_2, var = "X5.1482066461")
SNPclean3_x2 <-fisherpvalue(data1 = SNPclean3_1, data2 = SNPclean3_2, var = "X5.1482066462")

SNPBBclean3_1 <- mlr(names = names3, cov = SNPBBcov, data = miRNA_pheno3, var = "X5.1482066461:BBYes" )
SNPBBclean3_2 <- mlr(names = names3, cov = SNPBBcov, data = miRNA_pheno3_2, var = "X5.1482066461:BBYes" )
SNPBBclean3_x1 <-fisherpvalue(data1 = SNPBBclean3_1, data2 = SNPBBclean3_2, var = "X5.1482066461:BBYes")
SNPBBclean3_x2 <-fisherpvalue(data1 = SNPBBclean3_1, data2 = SNPBBclean3_2, var = "X5.1482066462:BBYes")

SNPclean3_1_p <- SNPclean3_x1 %>%  filter(.,p.value < 0.05)
SNPclean3_2_p <- SNPclean3_x2 %>%  filter(.,p.value < 0.05)
SNPBBclean3_1_p <- SNPBBclean3_x1 %>% filter(.,p.value < 0.05)
SNPBBclean3_2_p <- SNPBBclean3_x2 %>% filter(.,p.value < 0.05)


SNPclean_psig  <- unique(c(as.character(SNP_clean1_1_p$miRNA), as.character(SNP_clean2_1_p$miRNA), as.character(SNPclean3_1_p$miRNA), as.character(SNPclean3_2_p$miRNA)))
SNPBBclean_psig <- unique(c(as.character(SNPBB_clean1_1_p$miRNA), as.character(SNP_clean2_1_p$miRNA), as.character(SNPclean3_1_p$miRNA), as.character(SNPclean3_2_p$miRNA)))
SNPmlr_psig <- unique(c(as.character(SNPmlr1_1_p$miRNA), as.character(SNPmlr2_1_p$miRNA), as.character(SNPmlr3_1_p$miRNA), as.character(SNPmlr3_2_p$miRNA)))
SNPBBmlr_psig <- unique(c(as.character(SNPBBmlr1_1_p$miRNA), as.character(SNPBBmlr2_1_p$miRNA), as.character(SNPBBmlr3_1_p$miRNA), as.character(SNPBBmlr3_2_p$miRNA)))

SNPlist <- list(`SNP only` = SNPclean_psig,`SNP*BB only`=SNPBBclean_psig, `SNP + clinical`=SNPmlr_psig, `SNP*BB + clinical` = SNPBBmlr_psig )
venn.diagram(x = SNPlist, filename = "SNPVenndiagram.png",main = "SNP MLR Modeling", fill = c("red", "blue", "green", "purple"))

intersect(intersect(intersect(SNPclean_psig, SNPBBclean_psig), SNPmlr_psig), SNPBBmlr_psig)
SNPlist

summary(lm(miR_19a_3p~X10.115805056 + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno1))
summary(lm(miR_19a_3p~X10.115805056*BB + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno1))
summary(lm(miR_19a_3p~X10.115805056, data = miRNA_pheno1))
summary(lm(miR_19a_3p~X10.115805056*BB, data = miRNA_pheno1))
summary(lm(miR_186_5p_a2~X10.115805056 + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno1))
summary(lm(miR_186_5p_a2~X10.115805056*BB + AGE8 + SEX + HGT8 + WGT8 + rankcon + rankqual + rank260, data = miRNA_pheno1))
summary(lm(miR_186_5p_a2~X10.115805056, data = miRNA_pheno1))
summary(lm(miR_186_5p_a2~X10.115805056*BB, data = miRNA_pheno1))
