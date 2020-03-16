rm(list = ls())
library(emmeans)
library(tidyverse)
library(ggpubr)

setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")

topsnps<-read.csv("TopSNPsfinallist.csv", stringsAsFactors = F)
load("allSNPs_cleaned_andQC_genomewide.RData")

toppos <- topsnps$SNP
toppos <- as.integer(gsub("^..?:", "", toppos))
alldata <- tidy_clean2 %>% unnest(cols = c(data))
ADRB2data <- alldata %>% filter(., POS %in% toppos[topsnps$Gene == "ADRB2"]) 
ADRB2data_f <- ADRB2data[ADRB2data$SEX == "Female",]


ADRB2data_f$s8cbl24bd
ADRB2_Genotype_f <- lm(s8cbl24bd~Genotype*BB + AGE8 + HGT8 + EST8 + BMI8, data = ADRB2data_f)

p<- emmip(ADRB2_Genotype_f, formula = Genotype~BB) 
data <- p$data

ggplot(data = data, aes(x = BB, y = yvar, group = Genotype, color = Genotype)) + geom_point(position = position_dodge(0.9))+
  geom_errorbar(aes(ymin=yvar - SE, ymax=yvar + SE), width=.2, position=position_dodge(0.9)) + ylab("EMMean LS BMD") + xlab("BB user") + ggtitle("ADRB2 rs2400706 in Female Only")

emmip(ADRB2_Genotype_f, formula = Genotype~BB) + ylab("EMMean LS BMD") + xlab("BB user") + ggtitle("ADRB2 rs2400706 in Female Only")

table(ADRB2data$Genotype)
table(ADRB2data_f$Genotype)
mean(ADRB2data$AGE8[ADRB2data_f$Genotype == "normal"])
mean(ADRB2data$AGE8[ADRB2data_f$Genotype == "heterozygous"])
mean(ADRB2data$AGE8[ADRB2data_f$Genotype == "homozygous alternative"])

mean(ADRB2data_f$AGE8[ADRB2data_f$Genotype == "normal"])
mean(ADRB2data_f$AGE8[ADRB2data_f$Genotype == "heterozygous"])
mean(ADRB2data_f$AGE8[ADRB2data_f$Genotype == "homozygous alternative"])

table(ADRB2data_f$HRX8, ADRB2data_f$Genotype)
table(ADRB2data_f$HRX8, ADRB2data_f$Genotype)
