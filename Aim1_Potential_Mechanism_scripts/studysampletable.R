# Table 1: Characteristic of Study sample #
rm(list = ls())
# Install packages if you have not previously installed:
# install.packages("tidyverse")
#devtools::install_github("rstudio/gt")
#install.packages("remotes")
#remotes::install_github("ddsjoberg/gtsummary")

library(tidyverse)
library(gt)
library(dplyr)
library(gtsummary)
library(ggplot2)
library(emmeans)
library(reshape2)
library(matrixStats)
library(haven)
# Read in data
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
miRNAdat <- read.csv('l_mrna_2011_m_0797s_17_c1.csv')
pheno <- read.csv('PhenoData_5_28.csv')
mir_tech = "mirna_tech_17.sas7bdat"
mirnatech <- read_sas(mir_tech)

miRNA_delta_cq <- miRNAdat[-1]
miRNA_delta_cq <- -(miRNA_delta_cq-27)
miRNAdat <- cbind(miRNAdat[1], miRNA_delta_cq)

quantcon<-quantile(mirnatech$concentration, probs = seq(0,1, 0.1))
mirnatech$rankcon <- cut(mirnatech$concentration, quantcon)
quantqual <-quantile(mirnatech$RNA_quality, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rankqual <- cut(mirnatech$RNA_quality, quantqual)
quant260 <- quantile(mirnatech$`_260_280`, probs = seq(0,1, 0.1), na.rm = T)
mirnatech$rank260 <-cut(mirnatech$`_260_280`, quant260)
miRNA_pheno4 <- merge(pheno, mirnatech, by.x = "shareid", by.y = "shareid")
miRNA_pheno <- merge(miRNA_pheno4, miRNAdat, by.x = "shareid", by.y = "shareid")
miRNA_pheno$Isolation_Batch<-as.factor(miRNA_pheno$Isolation_Batch)



pheno1 <- pheno
pheno_female <- pheno1 %>% filter(., SEX == "Female")
pheno_male <- pheno1 %>%  filter(., SEX == "Male")
pheno1 %>% tbl_summary(.,by = "BB",statistic = list(..continuous.. = "{mean} ({sd})",..categorical.. = "{n} / {N} ({p}%)")) 
pheno1 %>% tbl_summary(., statistic = list(..continuous.. = "{mean} ({sd})",..categorical.. = "{n} / {N} ({p}%)")) 
pheno1$Tscore <- ((-0.023 + 0.939 * pheno1$f8cbnbmd  - 0.019)/1.087 - 0.858) / 0.120
t.test(f8cbnbmd~BB, pheno1)# Not significant
t.test(Tscore~BB, pheno1)# Not significant

t.test(f8cbtobmd~BB, pheno1)# Sig (p = 0.01)
t.test(f8cbtrbmd~BB, pheno1)# Sig (p = 0.004)
t.test(s8cbl2bd~BB, pheno1)# Sig (p = 0.0002)
t.test(s8cbl3bd~BB, pheno1)# Highly sig p < 0.0001
t.test(s8cbl4bd~BB, pheno1)# Highly sig p < 0.0001
t.test(s8cbl24bd~BB, pheno1)# Highly sig (p < 0.0001)

t.test(f8cbnbmd~BB, pheno_female)# Not significant
t.test(f8cbtobmd~BB, pheno_female)# Not significant
t.test(s8cbl24bd~BB, pheno_female)# Sig (p = 0.007)


t.test(f8cbnbmd~BB, pheno_male)# Not significant
t.test(f8cbtobmd~BB, pheno_male)# Not significant
t.test(s8cbl24bd~BB, pheno_male)# Sig (p = 0.01)

FToMeanBB <- mean(pheno1$f8cbtobmd[pheno1$BB == "Yes"], na.rm = T)
FToMeanBBNo <- mean(pheno1$f8cbtobmd[pheno1$BB == "No"], na.rm = T)
S24MeanBB <- mean(pheno1$s8cbl24bd[pheno1$BB == "Yes"], na.rm = T)
S24MeanBBNo <- mean(pheno1$s8cbl24bd[pheno1$BB == "No"], na.rm = T)

FToSEBB <- sd(pheno1$f8cbtobmd[pheno1$BB == "Yes"])/ sqrt(length(pheno1$f8cbtobmd[pheno1$BB == "Yes"]))
FToSEBBNo <- sd(pheno1$f8cbtobmd[pheno1$BB == "No"])/ sqrt(length(pheno1$f8cbtobmd[pheno1$BB == "No"]))
S24SEBB <- sd(pheno1$s8cbl24bd[pheno1$BB == "Yes"], na.rm = T)/ sqrt(length(pheno1$s8cbl24bd[pheno1$BB == "Yes"]))
S24SEBBNo <- sd(pheno1$s8cbl24bd[pheno1$BB == "No"], na.rm = T)/ sqrt(length(pheno1$s8cbl24bd[pheno1$BB == "No"]))

fnmodel <- lm(f8cbnbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
ftomodel<-lm(f8cbtobmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
ftrmodel <- lm(f8cbtrbmd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s2model <- lm(s8cbl2bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s3model <- lm(s8cbl3bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s4model <- lm(s8cbl4bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)
s24model <- lm(s8cbl24bd~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno1)

summary(fnmodel)
summary(ftomodel)
summary(ftrmodel)
summary(s2model)
summary(s3model)
summary(s4model)
summary(s24model)

emmeans(fnmodel, specs = "BB")
emmeans(ftomodel, specs = "BB")
fto<-emmeans(ftomodel, specs = "BB")
emmeans(ftrmodel, specs = "BB")
emmeans(s2model, specs = "BB")
emmeans(s3model, specs = "BB")
emmeans(s4model, specs = "BB")
s24<-emmeans(s24model, specs = "BB")
emmeans(s24model, specs = "BB")

fto <- data.frame(fto, BMD = "Total Femur")
s24 <- data.frame(s24, BMD = "Spine L2-L4")
bmd <- rbind(fto,s24)

FToMeasurements <- data.frame(BMD = pheno1$f8cbtobmd, BB = pheno1$BB, Group = "Total Femur")
LSMeasurements <- data.frame(BMD = pheno1$s8cbl24bd, BB = pheno1$BB, Group = "Spine LS-L4")
BMDMeasurements <- rbind(FToMeasurements, LSMeasurements)

write.csv(bmd, "BMDEMMeans.csv")
ggplot(data = bmd, aes(x=BMD, y=emmean, fill = BB))+ geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " BMD by BB use", x = "BMD Measurement Type", y = "BMD", fill = "BB User Status") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette="Paired")

ggplot(data = BMDMeasurements, aes(x=Group, y=BMD, fill = BB))+ geom_boxplot(position=position_dodge()) + 
  labs(title = " BMD by BB use", x = "BMD Measurement Type", y = "BMD", fill = "BB User Status") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette="Paired")

ggplot(data = bmd, aes(x=BMD, y=emmean, color = BB))+ geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " BMD by BB use", x = "BMD Measurement Type", y = "BMD", fill = "BB User Status") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")


# ggplot(data = bmd, aes(x=BB, y=emmean, fill = BMD))+ geom_dotplot(binaxis = "y", stackdir='center',  position=position_dodge()) + 
#   geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width=.2,position=position_dodge(.9)) + 
#   labs(title = " BMD by BB use", x = "BB User Status", y = "BMD", fill = "BMD Type") + theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette="Paired")


#boxplot(f8cbtobmd~BB, data = pheno1, xlab = "BB User",ylab = "Total Femur BMD", main = "Total Femur BMD by BB use" )
#boxplot(s8cbl24bd~BB, data = pheno1, xlab = "BB User",ylab = "Spine L2-L4 BMD", main = "Spine L2-L4 BMD by BB use" )

#FTo <- data.frame(BB_status = c("Yes", "No"), BMD_mean = c(FToMeanBB, FToMeanBBNo), BMD_SE = c(FToSEBB, FToSEBBNo))
#ggplot(data = FTo, aes(x=BB_status, y=BMD_mean))+ geom_bar(stat = "identity", fill = "gray") + 
#  geom_errorbar(aes(ymin = BMD_mean - BMD_SE, ymax = BMD_mean + BMD_SE)) + 
#  labs(title = "Total Femur BMD by BB use", x = "BB User Status", y = "Total Femur BMD") + theme_classic() +
#  theme(plot.title = element_text(hjust = 0.5))
 
#S24 <- data.frame(BB_status = c("Yes", "No"), BMD_mean = c(S24MeanBB, S24MeanBBNo), BMD_SE = c(S24SEBB, S24SEBBNo))
#ggplot(data = S24, aes(x=BB_status, y=BMD_mean))+ geom_bar(stat = "identity", fill = "gray") + 
#  geom_errorbar(aes(ymin = BMD_mean - BMD_SE, ymax = BMD_mean + BMD_SE)) + 
#  labs(title = "Spine L2-L4 BMD by BB use", x = "BB User Status", y = "Spine L2-L4 BMD") + theme_classic() +
#  theme(plot.title = element_text(hjust = 0.5))

plot(f8cbtobmd~AGE8, data = pheno1, main = "Total Femur BMD by Age", xlab = "Age at Exam 8", ylab ="Total Femur BMD")
abline(lm(f8cbtobmd~AGE8, data = pheno1), col="red") 
summary(lm(f8cbtobmd~AGE8, data = pheno1))
plot(s8cbl24bd~AGE8, data = pheno1, main = "Spine L2-L4 BMD by Age", xlab = "Age at Exam 8", ylab ="Spine L2-L4 BMD")
abline(lm(s8cbl24bd~AGE8, data = pheno1), col="red")
summary(lm(s8cbl24bd~AGE8, data = pheno1))

BB_pheno <- miRNA_pheno[miRNA_pheno$BB == "Yes",]
drop <- c("SEX", "CURRSMK8", "DMRX8", "HRX8", "LIPRX8", "EST8", "BB", "B1", "menov", "priorcvd")
BB_pheno1 <- BB_pheno[,!names(BB_pheno) %in% drop]
BB_pheno1 <- BB_pheno1[-c(1:12)]
noBB_pheno <- miRNA_pheno[miRNA_pheno$BB == "No",]
noBB_pheno1 <- noBB_pheno[,!names(noBB_pheno) %in% drop]
noBB_pheno1 <- noBB_pheno1[-c(1:12)]
miRNA_pheno_ave <-miRNA_pheno[,!names(miRNA_pheno) %in% drop] 
# miRNA_pheno_ave <- colMeans(miRNA_pheno_ave, na.rm = T)
# BB_means<-data.frame(var = names(BB_pheno1), Means = colMeans(BB_pheno1, na.rm = T), SE = colSds(as.matrix(BB_pheno1), na.rm = T)/sqrt(nrow(BB_pheno1)), group = "BB user")
# 
# noBB_means <- data.frame(var = names(noBB_pheno1), Means = colMeans(noBB_pheno1, na.rm = T), SE = colSds(as.matrix(BB_pheno1), na.rm = T)/sqrt(nrow(noBB_pheno1)), group = "non-BB user")
# Means <- rbind(BB_means, noBB_means)
library(pheatmap)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Clustering_miRNA/Figures")
WGCNA <- read.csv("MiRNA_significance_module_membership_nofilter.csv")
WGCNA_blue <- WGCNA %>% filter(.,mergedColors == "blue")
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
int <- read.csv("Spineoverlapresults_7_9_rank.csv")

# Means$var <- as.character(Means$var)
# Means_blue <- Means %>%  filter(., var %in% as.character(WGCNA_blue$X))
# Means_int <- Means %>% filter(., var %in% as.character(int$miRNA))
# pheatmap(na.omit(Means_blue[2:3]))
# pheatmap(na.omit(Means_blue[2:3]), scale = "row")
# pheatmap(na.omit(Means_int[2:3]))
# pheatmap(na.omit(Means_int[2:3]), scale = "row")

# write.csv(Means_int, "FigureS2data.csv")
# ggplot(data = Means_int, aes(x=var, y = Means, color = group)) + geom_point(position=position_dodge(width = 0.9)) + geom_errorbar(aes(ymin = Means - SE, ymax = Means + SE), width=.2, position = position_dodge(0.9))+theme_minimal()+scale_color_brewer(palette = "Paired")
boxplot(miR_19a_3p~BB,data = miRNA_pheno)
miRNA_pheno$miR19a <- miRNA_pheno$miR_19a_3p
miRNA_pheno$miR186 <- miRNA_pheno$miR_186_5p_a2
boxplot(miR19a~BB, data = miRNA_pheno, xlab = "BB User Status", ylab = "| \u0394 Cq |", main = "miR-19a-3p by BB Use")

model19 <- lm(miR19a~ f8cbtobmd + BB + f8cbtobmd*BB , data = miRNA_pheno) 
summary(model19)
model186 <- lm(miR186~ f8cbtobmd + BB + f8cbtobmd*BB , data = miRNA_pheno) 
summary(model186)
model186 <- lm(miR186~ BB, data = miRNA_pheno) 

ggplot(data = miRNA_pheno, aes(x = f8cbtobmd, y = miR19a, color = BB)) + geom_smooth(method = lm, se = TRUE, fullrange = TRUE) + ggtitle(label =  "miR-19a-3p by Total Femur BMD") + 
  xlab("Total Femur BMD") + ylab ("| \u0394 Cq |") + theme_minimal() + scale_color_brewer(palette = "Reds") + theme(plot.title = element_text(hjust = 0.5)) +theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10), legend.position = c(0.1,0.8), legend.background = element_rect(fill="white",colour = "white")) 
ggsave(filename = "miR19byTFbmd.tiff", width = 90, height = 60 ,unit = "mm")

ggplot(data = miRNA_pheno, aes(x = f8cbtobmd, y = miR186, color = BB)) + geom_smooth(method = lm, se = TRUE, fullrange = TRUE) + ggtitle(label =  "miR-186-5p by Total Femur BMD") + 
  xlab("Total Femur BMD") + ylab ("| \u0394 Cq |") + theme_minimal() + scale_color_brewer(palette = "Reds") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), title = element_text(size = 10), legend.text = element_text(size = 10), legend.position = c(0.1,0.8), legend.background = element_rect(fill="white",colour = "white")) 
ggsave(filename = "miR186byTFbmd.tiff", width = 90, height = 60 ,unit = "mm")

miRNA_pheno$FtoGroup<-ifelse(miRNA_pheno$f8cbtobmd<mean(miRNA_pheno$f8cbtobmd), "Low", "High")
model19a <- lm(miR19a ~ FtoGroup + BB + FtoGroup*BB, data = miRNA_pheno)
summary(model19a)
EM_19a <- as.data.frame(emmeans(model19a, pairwise ~ FtoGroup | BB))

ggplot(data =EM_19a, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "miR-19a-3p expression in BMD by BB use", x = "BB User Status", y = "EMMean of | \u0394 Cq |", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

model19aBB <- lm(miR19a ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(model19aBB)
EM_19a <- as.data.frame(emmeans(model19aBB, pairwise ~  BB))
summary(lm(miR19a~BB, data = miRNA_pheno))

model186BB <- lm(miR186 ~ BB + AGE8 + SEX + HGT8 + WGT8 + rank260 + rankcon + rankqual, data = miRNA_pheno)
summary(model186BB)
EM_186 <- as.data.frame(emmeans(model186BB, pairwise ~  BB))
