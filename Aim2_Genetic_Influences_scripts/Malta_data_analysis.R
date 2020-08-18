rm(list = ls())
setwd("~/")
library(xlsx)
options(stringsAsFactors = F)
data<-read.xlsx("Validation_Study_MOFS_July2020.xlsx", sheetIndex = 1)
names(data)
data <- data[-c(1059:1061),]

data$FN_BMD <- as.numeric(data$sBMD)
data$BB <- as.factor(data$BetaBlockers)
data$rs12414657..T.C. <- as.numeric(data$rs12414657..T.C.)
data$rs11124190..C.G. <- as.numeric(data$rs11124190..C.G.)
data$rs12414657..T.C.[data$rs12414657..T.C. == 0]<- NA

table(data$BB)
hist(data$FN_BMD)
hist(data$rs11124190..C.G.)
hist(data$rs12414657..T.C.)
hist(data$Age) #Not Normally distributed
hist(data$Height)
hist(data$BMI) #Long right tail

rs124_ADRB1<-lm("FN_BMD ~ BB*rs12414657..T.C. + Age + Height + BMI", data = data) # No estrogen use available
summary(rs124_ADRB1)
rs111_HDAC4 <- lm("FN_BMD ~ BB*rs11124190..C.G. + Age + Height + BMI", data = data) # No estrogen use available
summary(rs111_HDAC4)

