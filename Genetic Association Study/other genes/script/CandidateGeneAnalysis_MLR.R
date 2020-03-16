# Run Analysis
rm(list = ls())
library(tidyverse)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data/")
load("allSNPs_cleaned_andQC_genomewide.RData")
all_snps <- tidy_clean2

#Female Only and Male Only data
all_snps_F <- all_snps %>% unnest(cols = c(data)) %>% filter(., SEX == "Female") %>% group_by(., POS) %>% nest()
all_snps_M <- all_snps %>% unnest(cols = c(data)) %>% filter(., SEX == "Male") %>% group_by(., POS) %>% nest()

# MLR Function for SNPs
snp_mlr<-function(BMD, cov, data, length){
  mlr1 <- NULL
  for (i in 1:length){
    mlr <- lm(as.formula(paste(BMD, cov)), data = data[[2]][[i]])
    x <- summary(mlr)
    z <- as.data.frame(x$coefficients)
    z$SNP <- paste(unique(data[[2]][[i]]$Gene),data$POS[i], sep = ":")
    z$Gene <- unique(data[[2]][[i]]$Gene)
    z$var <- row.names(z)
    row.names(z)<- NULL
    if (i == 1){
      mlr1 <- z
    }else{
      mlr1 <- rbind(mlr1, z)
    }
  }
  return(mlr1)
}

snp_mlr_model<-function(BMD, cov, data, length){
  mlr1 <- list()
  for (i in 1:length){
    mlr <- lm(as.formula(paste(BMD, cov)), data = data[[2]][[i]])
    mlr1[[i]]<-mlr
  }
  return(mlr1)
}
# Sex Combined Model ####
# Total Femur (None significant)
To_F<- snp_mlr(BMD = "f8cbtobmd",cov = "~ Genotype*BB + AGE8 + SEX + HGT8 + EST8 + BMI8",data = all_snps, length = length(all_snps$data))
To_F$Model <- "ToF:SC"
sig_To_F_Bonferroni <- To_F %>% filter(.,`Pr(>|t|)` < 1.22971e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F <- To_F %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F_gene <- To_F %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F_05 <- To_F %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_To_F_05$SNP))
length(unique(sig_To_F$SNP))
length(unique(sig_To_F_gene$SNP))

# Femoral Neck (None significant)
FN_F<- snp_mlr(BMD = "f8cbnbmd",cov = "~ Genotype*BB + AGE8 + SEX + HGT8 + EST8 + BMI8",data = all_snps, length = length(all_snps$data))
FN_F$Model <- "FN:SC"
sig_FN_F_Bonferroni <- FN_F %>% filter(.,`Pr(>|t|)` < 1.22971e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_F <- FN_F %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_F_gene <- FN_F %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_F_gene4 <- FN_F %>% filter(.,`Pr(>|t|)` < 0.0008333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_F_05 <- FN_F %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_FN_F_05$SNP))
length(unique(sig_FN_F_05$Gene))
length(unique(sig_FN_F$SNP))
length(unique(sig_FN_F$Gene))
length(unique(sig_FN_F_gene$SNP))
length(unique(sig_FN_F_gene$Gene))
length(unique(sig_FN_F_gene4$SNP))
length(unique(sig_FN_F_gene4$Gene))

# Femoral Trochanter (None significant)
Tr_F<- snp_mlr(BMD = "f8cbtrbmd",cov = "~ Genotype*BB + AGE8 + SEX + HGT8 + EST8 + BMI8",data = all_snps, length = length(all_snps$data))
Tr_F$Model <- "Tr:SC"
sig_Tr_F_Bonferroni <- Tr_F %>% filter(.,`Pr(>|t|)` < 1.22971e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_Tr_F <- Tr_F %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_Tr_F_gene <- Tr_F %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_Tr_F_05 <- Tr_F %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_Tr_F_05$SNP))
length(unique(sig_Tr_F$SNP))
length(unique(sig_Tr_F_gene$SNP))

# Spine (None significant)
LS <- snp_mlr(BMD = "s8cbl24bd",cov = "~ Genotype + BB + Genotype*BB + AGE8 + SEX + HGT8 + EST8 + BMI8",data = all_snps, length = length(all_snps$data))
LS$Model <- "LS:SC"
sig_LS_Bonferroni <- LS %>% filter(.,`Pr(>|t|)` < 1.22971e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS <- LS %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_gene <- LS %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_gene4 <- LS %>% filter(.,`Pr(>|t|)` < 0.0008333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_05 <- LS %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_LS_05$SNP))
length(unique(sig_LS_05$Gene))
length(unique(sig_LS$SNP))
length(unique(sig_LS$Gene))
length(unique(sig_LS_gene$SNP))
length(unique(sig_LS_gene$Gene))
length(unique(sig_LS_gene4$SNP))
length(unique(sig_LS_gene4$Gene))

#Female Only Model ####
# Total Femur (None significant)
To_F_Fem<- snp_mlr(BMD = "f8cbtobmd",cov = "~ Genotype*BB + AGE8 + HGT8 + EST8 + BMI8",data = all_snps_F, length = length(all_snps_F$data))
To_F_Fem$Model <- "ToF:F"
sig_To_F_Fem_Bonferroni <- To_F_Fem %>% filter(.,`Pr(>|t|)` < 1.226091e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F_Fem <- To_F_Fem %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F_Fem_gene <- To_F_Fem %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_To_F_05_Fem <- To_F_Fem %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_To_F_05_Fem$SNP))
length(unique(sig_To_F_Fem$SNP))
length(unique(sig_To_F_Fem_gene$SNP))

# Femoral Neck (None significant)
FN_Fem<- snp_mlr(BMD = "f8cbnbmd",cov = "~ Genotype*BB + AGE8 + HGT8 + EST8 + BMI8",data = all_snps_F, length = length(all_snps_F$data))
FN_Fem$Model <- "FN:F"
sig_FN_Fem_Bonferroni <- FN_Fem %>% filter(.,`Pr(>|t|)` < 1.226091e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_Fem <- FN_Fem %>% filter(.,`Pr(>|t|)` < 5.170631e-05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_Fem_gene <- FN_Fem %>% filter(.,`Pr(>|t|)` < 0.003333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_Fem_gene4 <- FN_Fem %>% filter(.,`Pr(>|t|)` < 0.0008333333) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FN_Fem_05 <- FN_Fem %>% filter(.,`Pr(>|t|)` < 0.05) %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_FN_Fem_05$SNP))
length(unique(sig_FN_Fem_05$Gene))
length(unique(sig_FN_Fem$SNP))
length(unique(sig_FN_Fem$Gene))
length(unique(sig_FN_Fem_gene$SNP))
length(unique(sig_FN_Fem_gene$Gene))
length(unique(sig_FN_Fem_gene4$SNP))
length(unique(sig_FN_Fem_gene4$Gene))

# Femoral Trochanter (None significant)
FTr_Fem<- snp_mlr(BMD = "f8cbtrbmd",cov = "~ Genotype + BB + Genotype*BB + AGE8 + HGT8 + EST8 + BMI8",data = all_snps_F, length = length(all_snps_F$data))
FTr_Fem$Model <- "Tr:F"
sig_FTr_Fem_Bonferroni <- FTr_Fem %>% filter(.,`Pr(>|t|)` < 1.226091e-05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FTr_Fem <- FTr_Fem %>% filter(.,`Pr(>|t|)` < 5.170631e-05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FTr_Fem_gene <- FTr_Fem %>% filter(.,`Pr(>|t|)` < 0.003333333)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_FTr_05_Fem <- FTr_Fem %>% filter(.,`Pr(>|t|)` < 0.05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_FTr_05_Fem$SNP))
length(unique(sig_FTr_Fem$SNP))
length(unique(sig_FTr_Fem_gene$SNP))

# Spine (None significant)
LS_Fem <- snp_mlr(BMD = "s8cbl24bd",cov = "~ Genotype + BB + Genotype*BB + AGE8 + HGT8 + EST8 + BMI8",data = all_snps_F, length = length(all_snps_F$data))
LS_Fem$Model <- "LS:F"
sig_LS_Fem_Bonferroni <- LS_Fem %>% filter(.,`Pr(>|t|)` < 1.226091e-05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_Fem <- LS_Fem %>% filter(.,`Pr(>|t|)` < 5.170631e-05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_Fem_gene <- LS_Fem %>% filter(.,`Pr(>|t|)` < 0.003333333)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_Fem_gene4 <- LS_Fem %>% filter(.,`Pr(>|t|)` < 0.0008333333)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
sig_LS_Fem_05 <- LS_Fem %>% filter(.,`Pr(>|t|)` < 0.05)  %>% filter(., var %in% c("Genotypeheterozygous:BBYes", "Genotypehomozygous alternative:BBYes"))
length(unique(sig_LS_Fem_05$SNP))
length(unique(sig_LS_Fem_05$Gene))
length(unique(sig_LS_Fem$SNP))
length(unique(sig_LS_Fem$Gene))
length(unique(sig_LS_Fem_gene$SNP))
length(unique(sig_LS_Fem_gene$Gene))
length(unique(sig_LS_Fem_gene4$SNP))
length(unique(sig_LS_Fem_gene4$Gene))

simplem_results <- rbind(sig_FN_F, sig_FN_Fem, sig_FTr_Fem, sig_LS, sig_LS_Fem, sig_To_F, sig_To_F_Fem, sig_Tr_F)
gene_results <- rbind(sig_FN_F_gene, sig_FN_Fem_gene, sig_FTr_Fem_gene, sig_LS_Fem_gene, sig_LS_gene, sig_To_F_Fem_gene, sig_To_F_Fem_gene, sig_Tr_F_gene)
p_results <- rbind(sig_FN_05_Fem, sig_FN_F_05, sig_FTr_05_Fem, sig_LS_05, sig_LS_05_Fem, sig_To_F_05, sig_To_F_05_Fem, sig_Tr_F_05)
all_results <- rbind(FN_F, FN_Fem, To_F, To_F_Fem, Tr_F, FTr_Fem, LS, LS_Fem)

write.csv(simplem_results, "SigSNPsResults_simpleM.csv", quote = F, row.names = F)
write.csv(gene_results, "SigSNPsResults_gene.csv", quote = F, row.names = F)
write.csv(p_results, "SigSNPsResults_p05.csv", quote = F, row.names = F)
write.csv(all_results, "AllMLR_Results.csv", quote = F, row.names = F)

table(simplem_results$Gene)
length(unique(simplem_results$SNP))
table(gene_results$Gene)
length(unique(gene_results$SNP))
table(p_results$Gene)
length(unique(p_results$SNP))
