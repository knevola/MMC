#ESR1, MAPK3, EGF,IGF1, STK11, IGFBP1, PSTPIP1, MECP2, WNK1, ZNF446
rm(list =ls())
library(readr)
library(emmeans)
library(ggplot2)
setwd('/home/clary@mmcf.mehealth.org/Framingham/OmicData/data')
GPL5175 <- read_table2("GPL5175.txt", skip = 14)
mRNA_sample <- read.delim("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data/phe000002.v7_release_manifest.txt", comment.char="#")
mRNAdat <- read_delim("FinalFile_Gene_OFF_2446_Adjusted_c1.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pheno <- read.csv('PhenoData_5_28.csv')
pheno_samp <- merge(pheno, mRNA_sample[,c(1,2)], by.x = "shareid", by.y = "Subject_ID")
data_adj_t <-as.data.frame(t(mRNAdat))
data_adj_t <- data_adj_t[-1,]
colnames(data_adj_t)<- mRNAdat$transcript_cluster_id
dat<- merge(pheno_samp, data_adj_t, by.x = "Sample_ID", by.y = 0)
dat$WNK1 <- dat$`3400034`
dat$ESR1 <- dat$`2931763`
dat$MAPK3 <- dat$`3687494`
dat$EGF <- dat$`2739308`
dat$IGF1 <- dat$`3468345`
dat$STK11 <- dat$`3815566`
dat$IGFBP1 <- dat$`3000503`
dat$PSTPIP1 <- dat$`3602767`
dat$MECP2 <- dat$`4027056`
dat$ZNF446 <- dat$`3844175`
dat$IL1B <- dat$`2571510`
dat$RANKL <- dat$`3487299`
dat$ADRB1 <- dat$`3265140`
dat$ADRB2 <- dat$`2834743`

dat$FtoGroup<-ifelse(dat$f8cbtobmd>mean(dat$f8cbtobmd), "High", "Low")
dat$LSGroup <-ifelse(dat$s8cbl24bd>mean(dat$s8cbl24bd, na.rm = T), "High", "Low")

modelWNK1 <- aov(WNK1 ~ FtoGroup + BB + FtoGroup*BB, data = dat)
summary(modelWNK1)
meansWNK1 <- as.data.frame(emmeans(modelWNK1, pairwise ~ FtoGroup | BB))

ggplot(data = meansWNK1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " WNK1 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")
modelWNK1_clinical <- lm(WNK1 ~ FtoGroup + BB + FtoGroup*BB + AGE8 + SEX + HGT8 + WGT8, data = dat)
summary(modelWNK1_clinical)
meansWNK1_clinical <- as.data.frame(emmeans(modelWNK1_clinical, pairwise ~ FtoGroup | BB))

ggplot(data = meansWNK1_clinical, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " WNK1 expression in BMD by BB use with Clinical", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

#ESR1, MAPK3, EGF,IGF1, STK11, IGFBP1, PSTPIP1, MECP2, WNK1, ZNF446
modelESR1 <- aov(ESR1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelESR1) # NS
EM_ESR1<-as.data.frame(emmeans(modelESR1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_ESR1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " ESR1 expression in BMD by BB use (NS)", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelMAPK3 <- aov(MAPK3 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelMAPK3) # NS
EM_MAPK3<-as.data.frame(emmeans(modelMAPK3, pairwise ~ FtoGroup | BB))

ggplot(data =EM_MAPK3, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " MAPK3 expression in BMD by BB use (NS)", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")


modelEGF <- aov(EGF ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelEGF) # barely significant
EM_EGF <- as.data.frame(emmeans(modelEGF, pairwise ~ FtoGroup | BB))

ggplot(data =EM_EGF, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " EGF expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelIGF1 <- aov(IGF1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelIGF1) # BB
EM_IGF1 <- as.data.frame(emmeans(modelIGF1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_IGF1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = " IGF1 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelSTK11 <- aov(STK11 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelSTK11) # Fto
EM_STK11 <- as.data.frame(emmeans(modelSTK11, pairwise ~ FtoGroup | BB))

ggplot(data =EM_STK11, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "STK11 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelIGFBP1 <- aov(IGFBP1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelIGFBP1) # Fto
EM_IGFBP1 <- as.data.frame(emmeans(modelIGFBP1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_IGFBP1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "IGFBP1 expression in BMD by BB use (NS)", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelPSTPIP1 <- aov(PSTPIP1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelPSTPIP1) # Fto
EM_PSTPIP1 <- as.data.frame(emmeans(modelPSTPIP1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_PSTPIP1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "PSTPIP1 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelMECP2 <- aov(MECP2 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelMECP2) # Fto
EM_MECP2 <- as.data.frame(emmeans(modelMECP2, pairwise ~ FtoGroup | BB))

ggplot(data =EM_MECP2, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "MECP2 expression in BMD by BB use (NS)", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelWNK1 <- aov(WNK1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelWNK1) # Fto
EM_WNK1 <- as.data.frame(emmeans(modelWNK1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_WNK1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "WNK1 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")


modelZNF446 <- aov(ZNF446 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelZNF446) # Fto
EM_ZNF446 <- as.data.frame(emmeans(modelZNF446, pairwise ~ FtoGroup | BB))

ggplot(data =EM_ZNF446, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "ZNF446 expression in BMD by BB use(NS)", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelIL1B <- aov(IL1B ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelIL1B) # Fto
EM_IL1B <- as.data.frame(emmeans(modelIL1B, pairwise ~ FtoGroup | BB))

ggplot(data =EM_IL1B, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "IL1B expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelRANKL <- aov(RANKL ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelRANKL) # Fto
EM_RANKL <- as.data.frame(emmeans(modelRANKL, pairwise ~ FtoGroup | BB))

ggplot(data =EM_RANKL, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "RANKL expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelADRB1 <- aov(ADRB1 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelADRB1) # Fto
EM_ADRB1 <- as.data.frame(emmeans(modelADRB1, pairwise ~ FtoGroup | BB))

ggplot(data =EM_ADRB1, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "ADRB1 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")

modelADRB2 <- aov(ADRB2 ~ FtoGroup + BB + FtoGroup*BB , data = dat)
summary(modelADRB2) # Fto
EM_ADRB2 <- as.data.frame(emmeans(modelADRB2, pairwise ~ FtoGroup | BB))

ggplot(data =EM_ADRB2, aes(x = emmeans.BB, y = emmeans.emmean, color = emmeans.FtoGroup)) + geom_point(position=position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = emmeans.emmean - emmeans.SE, ymax = emmeans.emmean + emmeans.SE), width=.2,position=position_dodge(.9)) + 
  labs(title = "ADRB2 expression in BMD by BB use", x = "BB User Status", y = "EMMean", color = "FtoGroup") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_brewer(palette = "Paired")
