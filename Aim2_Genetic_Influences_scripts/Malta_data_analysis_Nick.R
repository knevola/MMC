rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")
library(xlsx)
options(stringsAsFactors = F)
data<-read.xlsx("HCS dataset (Pharmacogenetic Study of Beta Blocker Response in Bone).xlsx", sheetIndex = 2)
names(data)
data$FN_BMD <- as.numeric(data$bnnekbmd)
data$BB <- as.factor(data$hdrugc6)
data$est <- as.factor(data$hdrugh7)
data$rs12414657 <- as.numeric(data$rs12414657)
data$rs11124190 <- as.numeric(data$rs11124190)
data$rs34170507 <- as.numeric(data$rs34170507)
data$rs6567268 <- as.numeric(data$rs6567268)
data$rs12414657 <- 4 - data$rs12414657
data$rs34170507 <- 4 - data$rs34170507
fulldata <- data # hold full data set
# start with females
females <- data[data$absex==2,]
data <- females
table(data$BB)
hist(data$FN_BMD)
hist(data$rs11124190)
hist(data$rs12414657)
hist(data$bnage)
hist(data$aht)
hist(data$abmi)
table(data$est,exclude=F)

rs124_ADRB1<-lm("FN_BMD ~ BB*rs12414657 + bnage + aht + abmi + est", data = data)
summary(rs124_ADRB1)
rs111_HDAC4 <- lm("FN_BMD ~ BB*rs11124190 + bnage + aht + abmi + est", data = data)
summary(rs111_HDAC4)

# next males
males <- fulldata[fulldata$absex==1,]
data <- males
table(data$BB)
hist(data$FN_BMD)
hist(data$rs34170507)
hist(data$rs6567268)
hist(data$bnage)
hist(data$aht)
hist(data$abmi)

rs34170507_RANK <- lm("FN_BMD ~ BB*rs34170507 + bnage + aht + abmi", data = data)
summary(rs34170507_RANK)
rs6567268_RANK <- lm("FN_BMD ~ BB*rs6567268 + bnage + aht + abmi", data = data)
summary(rs6567268_RANK)

