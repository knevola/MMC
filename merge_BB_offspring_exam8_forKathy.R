rm(list=ls())
library(R.utils)
library(tidyverse)
library(tableone)
# offspring exam dates
datef = "phs000007.v30.pht003099.v5.p11.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt.gz"
datef = substring(datef,1,nchar(datef)-3)
date = read.delim(datef,skip=10,header=T,stringsAsFactors = F)
table(date$att8)
date8 = date %>% # 3795, 2803 who attended exam 8
  select(shareid,age8,att8,date8,idtype) %>% 
  filter(idtype==1) %>% 
  filter(!is.na(att8)) %>% 
  filter(att8==1)
# covariates
covf = "phs000007.v30.pht006027.v2.p11.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt.gz"
covf = substring(covf,1,nchar(covf)-3)
cov = read.delim(covf,skip=10,header=T,stringsAsFactors = F)
cov8 = cov %>% # 3795
  select(shareid,SEX,AGE8,BMI8,HGT8,WGT8,CURRSMK8,CPD8,DMRX8,HRX8,LIPRX8,IDTYPE) %>% 
  filter(IDTYPE==1)
estf = "phs000007.v30.pht000307.v8.p11.c1.meno1_8s.HMB-IRB-MDS.txt.gz"
estf = substring(estf,1,nchar(estf)-3)
est = read.delim(estf,skip=10,header=T,stringsAsFactors = F)
est8 = est %>%
  select(shareid,OVREM8,EST8,STOP_AGE)
# soe data
soef = "phs000007.v30.pht003335.v6.p11.c1.vr_soe4srv_2014_a_1027s.HMB-IRB-MDS.txt.gz"
soef = substring(soef,1,nchar(soef)-3)
soe = read.delim(soef,skip=10,header=T,stringsAsFactors = F) # 5184 in offspring
soe8 = soe[soe$IDTYPE==1,]
length(unique(soe$shareid[soe$IDTYPE==1])) # 1748
# do merge
cohort = merge(date8,cov8,by.x="shareid",by.y="shareid")
cohort = merge(cohort,est8,by.x="shareid",by.y="shareid",all.x=T)
cohort = cohort %>%
  mutate(menov = (OVREM8==2 | STOP_AGE<age8))
cohort$menov[cohort$SEX==1] = 2
cohort$EST8[cohort$SEX==1] = 2
cohortsoe = merge(cohort,soe,by.x="shareid",by.y="shareid",all.x=T)
getcvd = cohortsoe %>%
  mutate(cvd_event = ifelse(is.na(EVENT),0,ifelse(EVENT %in% c(1:19,30:49),1,0))) %>%
  mutate(event_days8 = DATE-date8) %>%
  mutate(priorcvd_event = ifelse(event_days8<0 & cvd_event==1,1,0)) %>%
  group_by(shareid) %>%
  summarize(priorcvd = max(priorcvd_event))
cohort = merge(cohort,getcvd,by.x="shareid",by.y="shareid",all.x=T)
cohort$priorcvd[is.na(cohort$priorcvd)] = 0
getdeath = cohortsoe %>%
  select(shareid,EVENT,DATE,date8) %>%
  mutate(death = ifelse(is.na(EVENT),0,ifelse(EVENT %in% c(21:29),1,0))) %>%
  mutate(censor_days8 = DATE-date8) %>%
  filter(censor_days8>0) %>%
  filter(death==1)
cohort = merge(cohort,getdeath,by.x="shareid",by.y="shareid",all.x=T)
cohort$death[is.na(cohort$death)] = 0
# offspring drug data
drugf = "phs000007.v30.pht000828.v6.p11.c1.meds1_8s.HMB-IRB-MDS.txt.gz"
drugf = substring(drugf,1,nchar(drugf)-3)
drug = read.delim(drugf,skip=10,header=T,stringsAsFactors = F)
drug = drug %>% 
  select(shareid,ther_gp1,MEDNAME,MEDPRN,chem_gp1,chem_nm1) %>% 
  filter(ther_gp1=="BETA BLOCKING AGENTS") %>% 
  filter(!(MEDNAME %in% c("ISTALOL","IC ATENOLOL"))) %>% 
  filter(MEDPRN==0) %>% # 4 RN, 50 unknown
  mutate(B1drug=ifelse(chem_gp1=="Beta blocking agents, selective",1,0)) %>% 
  mutate(BBdrug=1)  %>% 
  group_by(shareid) %>%
  summarize(BB=max(BBdrug),B1=max(B1drug))
cohort = merge(cohort,drug,by.x="shareid",by.y="shareid",all.x=T)
cohort$BB[is.na(cohort$BB)] = 0
cohort$B1[is.na(cohort$B1)] = 0
cohort = cohort %>%
  mutate(BBcat = ifelse(BB==0,"none",ifelse(B1==1,"selective","non-selective")))
cohort$BBcat = factor(cohort$BBcat,ordered=T,levels=c("none","non-selective","selective"))
# offspring BMD data
bmdofff = "phs000007.v30.pht003096.v2.p11.c1.t_bmdhs_2008_1_0748s.HMB-IRB-MDS.txt.gz"
bmdofff = substring(bmdofff,1,nchar(bmdofff)-3)
bmdoff = read.delim(bmdofff,skip=10,header=T,stringsAsFactors = F)
bmdoff = merge(bmdoff,date8,by.x="shareid",by.y="shareid")
bmd = bmdoff %>%
  mutate(f_exam8 = f8cbscdt-date8) %>%
  mutate(s_exam8 = s8cbscdt-date8) %>%
  mutate(mis_f8 = is.na(f8cbscdt) | f_exam8<0) %>%
  mutate(mis_s8 = is.na(s8cbscdt) | s_exam8<0) %>%
  filter(!mis_f8 | !mis_s8)
cohortbmd = merge(cohort,bmd,by.x="shareid",by.y="shareid")
# PS model for BB use
# you repeat for being in case control study
#modBB = glm(BB~AGE8+BMI8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmd,
            #family="binomial",na.action=na.exclude)
modBB = glm(BB~BMI8,data=cohortbmd,
            family="binomial",na.action=na.exclude)
summary(modBB)
score = predict(modBB,type="response")
unweighted = lm(f8cbnbmd~BB,data=cohortbmd)
summary(unadjusted)
notmissing = which(!is.na(score))
data2 = cohortbmd[notmissing,]
weights = 1/score[notmissing]
weighted = lm(f8cbnbmd~BB,data=data2,weights=weights)
summary(weighted)
adjusted = lm(f8cbnbmd~BB+AGE8+BMI8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmd)
summary(adjusted)
bmi_unweighted= aov(BMI8~factor(BB),data=data2)
summary(bmi_unweighted)
emmeans(bmi_unweighted,~BB)
bmi_weighted= aov(BMI8~factor(BB),data=data2,weights=weights)
summary(bmi_weighted)
emmeans(bmi_weighted,~BB)
hist(score[cohortbmd$BB==1])
hist(score[cohortbmd$BB==0],add=T,col="blue")

