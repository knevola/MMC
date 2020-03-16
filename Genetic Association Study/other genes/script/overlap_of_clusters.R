rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(tidyverse)
power5 <- read.csv("MiRNA_significance_module_membership_nofiltergenetics5.csv")
power8 <- read.csv("MiRNA_significance_module_membership_nofiltergenetics8.csv")
power10 <- read.csv("MiRNA_significance_module_membership_nofiltergenetics10.csv")
miRNA <- read.csv("MiRNA_significance_module_membership_nofilter_miRNA_only.csv")

names(power8) <- paste(names(power8), "_8", sep = "")
names(power10) <- paste(names(power10), "_10", sep = "")
names(power5) <- paste(names(power5), "_5", sep = "")
power8$mergedColors_8 <- paste(power8$mergedColors_8, "_8", sep = "")
power5$mergedColors_5 <- paste(power5$mergedColors_5, "_5", sep = "")
power10$mergedColors_10 <- paste(power10$mergedColors_10, "_10", sep = "")

miRNA_all <- merge(miRNA, power8, by.x = "miRNA", by.y = "miRNA_8")
miRNA_all <- merge(miRNA_all, power10, by.x = "miRNA", by.y = "miRNA_10")
miRNA_all <- merge(miRNA_all, power5, by.x = "miRNA", by.y = "miRNA_5")
miRNA_colors <- miRNA_all %>% select(., miRNA, mergedColors, mergedColors_5, mergedColors_8, mergedColors_10)

oby<- as.data.frame(table(miRNA_colors$mergedColors, miRNA_colors$mergedColors_5, miRNA_colors$mergedColors_8, miRNA_colors$mergedColors_10))
names(oby)<- c("Original", "Power_5", "Power_8", "Power_10", "Freq")


blue_miRNA <- miRNA_colors[miRNA_colors$mergedColors == "blue",]
blue_8<- as.data.frame(table(blue_miRNA$mergedColors_8))
total_8 <- as.data.frame(table(miRNA_colors$mergedColors_8))
new <- data.frame(setdiff(total_8$Var1, blue_8$Var1), 0)
names(new)<- names(blue_8)
blue_8 <- rbind(blue_8, new)
blue_8$Total <- total_8$Freq
blue_8$Percent_8 <- blue_8$Freq/blue_8$Total
blue_8$Percent_blue <- blue_8$Freq/43

blue_10<- as.data.frame(table(blue_miRNA$mergedColors_10))
total_10 <- as.data.frame(table(power10$mergedColors_10))
blue_10$Total <- total_10$Freq
blue_10$Percent_10 <- blue_10$Freq/blue_10$Total
blue_10$Percent_blue <- blue_10$Freq/43

blue_5<- as.data.frame(table(blue_miRNA$mergedColors_5))
total_5 <- as.data.frame(table(power5$mergedColors_5))
new <- data.frame(setdiff(total_5$Var1, blue_5$Var1), 0)
names(new)<- names(blue_5)
blue_5 <- rbind(blue_5, new)
blue_5$Total <- total_5$Freq
blue_5$Percent_5 <- blue_5$Freq/blue_5$Total
blue_5$Percent_blue <- blue_5$Freq/43

write.csv(blue_8, "Blue8.csv", row.names = F, quote = F)
write.csv(blue_10, "Blue10.csv", row.names = F, quote = F)
write.csv(blue_5, "Blue5.csv", row.names = F, quote = F)


brown_miRNA <- miRNA_colors[miRNA_colors$mergedColors == "brown",]
brown_8<- as.data.frame(table(brown_miRNA$mergedColors_8))
total_8 <- as.data.frame(table(miRNA_colors$mergedColors_8))
new <- data.frame(setdiff(total_8$Var1, brown_8$Var1), 0)
names(new)<- names(brown_8)
brown_8 <- rbind(brown_8, new)
brown_8$Total <- total_8$Freq
brown_8$Percent_8 <- brown_8$Freq/brown_8$Total
brown_8$Percent_brown <- brown_8$Freq/30

brown_10<- as.data.frame(table(brown_miRNA$mergedColors_10))
total_10 <- as.data.frame(table(power10$mergedColors_10))
new <- data.frame(setdiff(total_10$Var1, brown_10$Var1), 0)
names(new)<- names(brown_10)
brown_10 <- rbind(brown_10, new)
brown_10$Total <- total_10$Freq
brown_10$Percent_10 <- brown_10$Freq/brown_10$Total
brown_10$Percent_brown <- brown_10$Freq/30

brown_5<- as.data.frame(table(brown_miRNA$mergedColors_5))
total_5 <- as.data.frame(table(power5$mergedColors_5))
new <- data.frame(setdiff(total_5$Var1, brown_5$Var1), 0)
names(new)<- names(brown_5)
brown_5 <- rbind(brown_5, new)
brown_5$Total <- total_5$Freq
brown_5$Percent_5 <- brown_5$Freq/brown_5$Total
brown_5$Percent_brown <- brown_5$Freq/30

write.csv(brown_8, "brown8.csv", row.names = F, quote = F)
write.csv(brown_10, "brown10.csv", row.names = F, quote = F)
write.csv(brown_5, "brown5.csv", row.names = F, quote = F)

oby<- as.data.frame(table(miRNA_colors$mergedColors, miRNA_colors$mergedColors_5, miRNA_colors$mergedColors_8, miRNA_colors$mergedColors_10))
oby <- oby[oby$Freq != 0,]
oby5 <- as.data.frame(table(miRNA_colors$mergedColors, miRNA_colors$mergedColors_5))
oby5 <- oby5[oby5$Freq != 0,]
oby8 <- as.data.frame(table(miRNA_colors$mergedColors, miRNA_colors$mergedColors_8))
oby8 <- oby8[oby8$Freq != 0,]
oby10 <- as.data.frame(table(miRNA_colors$mergedColors, miRNA_colors$mergedColors_10))
oby10 <- oby10[oby10$Freq != 0,]

write.csv(oby, "ClusterChangeAcrossAllModels.csv", row.names = F, quote = F)
write.csv(oby5, "ClusterChangeAcrossOriginaland5.csv", row.names = F, quote = F)
write.csv(oby8, "ClusterChangeAcrossOriginaland8.csv", row.names = F, quote = F)
write.csv(oby10, "ClusterChangeAcrossOriginaland10.csv", row.names = F, quote = F)

oby5810<- as.data.frame(table(miRNA_colors$mergedColors_5, miRNA_colors$mergedColors_8, miRNA_colors$mergedColors_10))
oby5810 <- oby5810[oby5810$Freq != 0,]
oby58 <- as.data.frame(table(miRNA_colors$mergedColors_5, miRNA_colors$mergedColors_8))
oby58 <- oby58[oby58$Freq != 0,]
oby510 <- as.data.frame(table(miRNA_colors$mergedColors_5, miRNA_colors$mergedColors_10))
oby510 <- oby510[oby510$Freq != 0,]
oby810 <- as.data.frame(table(miRNA_colors$mergedColors_8, miRNA_colors$mergedColors_10))
oby810 <- oby810[oby810$Freq != 0,]

write.csv(oby5810, "ClusterChangeAcrossModels5810.csv", row.names = F, quote = F)
write.csv(oby58, "ClusterChangeAcrossModels58.csv", row.names = F, quote = F)
write.csv(oby510, "ClusterChangeAcrossModels510.csv", row.names = F, quote = F)
write.csv(oby810, "ClusterChangeAcrossModels810.csv", row.names = F, quote = F)

setwd("/home/clary@mmcf.mehealth.org/Framingham/eQTL Analysis/data")
miRNAs <- read.csv("SignificantMiRNAswithGenes.csv", stringsAsFactors = F)
RANK_miRNAs <- miRNAs[miRNAs$Var1 == "RANK",]
ADRB2_miRNAs <- miRNAs[miRNAs$Var1 == "ADRB2",]
FoxO1_miRNAs <- miRNAs[miRNAs$Var1 == "FoxO1",]

sig_miRNAs <- miRNA_colors %>% filter(., miRNA %in% miRNAs$Var2)
sig_miRNAs_RANK <- miRNA_colors %>% filter(., miRNA %in% RANK_miRNAs$Var2)
sig_miRNAs_ADRB2 <- miRNA_colors %>% filter(., miRNA %in% ADRB2_miRNAs$Var2)
sig_miRNAs_FoxO1 <- miRNA_colors %>% filter(., miRNA %in% FoxO1_miRNAs$Var2)

table(sig_miRNAs$mergedColors)
table(sig_miRNAs$mergedColors_5)
table(sig_miRNAs$mergedColors_8)
table(sig_miRNAs$mergedColors_10)

table(sig_miRNAs_RANK$mergedColors)
table(sig_miRNAs_RANK$mergedColors_5)
table(sig_miRNAs_RANK$mergedColors_8)
table(sig_miRNAs_RANK$mergedColors_10)

table(sig_miRNAs_ADRB2$mergedColors)
table(sig_miRNAs_ADRB2$mergedColors_5)
table(sig_miRNAs_ADRB2$mergedColors_8)
table(sig_miRNAs_ADRB2$mergedColors_10)

table(sig_miRNAs_FoxO1$mergedColors)
table(sig_miRNAs_FoxO1$mergedColors_5)
table(sig_miRNAs_FoxO1$mergedColors_8)
table(sig_miRNAs_FoxO1$mergedColors_10)
