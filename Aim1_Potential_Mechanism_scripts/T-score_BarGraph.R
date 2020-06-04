rm(list = ls())
library(tidyverse)
library(emmeans)
library(ggplot2)

setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data')
pheno <- read.csv('PhenoData_5_28.csv')

# Determine T-score differences between BB user vs non-users ####
pheno$FN_Tscore <- ((-0.023 + 0.939 * pheno$f8cbnbmd  - 0.019)/1.087 - 0.858) / 0.120

mean(pheno$FN_Tscore)
sd(pheno$FN_Tscore)

mean(pheno$FN_Tscore[pheno$BB == "Yes"])
sd(pheno$FN_Tscore[pheno$BB == "Yes"])

mean(pheno$FN_Tscore[pheno$BB == "No"])
sd(pheno$FN_Tscore[pheno$BB == "No"])

fnmodel <- lm(FN_Tscore~BB + AGE8 + SEX + HGT8 + WGT8, data = pheno)
summary(fnmodel)
fn<-as.data.frame(emmeans(fnmodel, specs = "BB"))

p<-ggplot(data=fn, aes(x = BB, y = emmean, fill = BB)) + geom_bar(stat = "identity") +theme_minimal()

p <- p + scale_fill_brewer(palette = "Blues")+ xlab("BB use") + ylab("Adjusted T-score")
p + geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2,
                  position=position_dodge(.9))
