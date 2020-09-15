rm(list=ls())
library(R.utils)
library(tidyverse)
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
# offspring fx data
fxofff = "phs000007.v30.pht001044.v7.p11.c1.vr_fxrev_2013_1_0847s.HMB-IRB-MDS.txt.gz"
fxofff = substring(fxofff,1,nchar(fxofff)-3)
fxoff = read.delim(fxofff,skip=10,header=T,stringsAsFactors = F)
fxoff = merge(fxoff,date8)
fx = fxoff %>% 
  select(shareid,of_FxSite,of_fxcircum,of_FxDate,of_fxside,date8,of_fxcircum) %>% 
  mutate(fx_days8 = of_FxDate-date8) %>%
  filter(fx_days8>0) %>%
  filter(!(of_FxSite %in% c(4,11,21))) %>%  #4: Skull or facial bones: includes jaw, nose, cheek, 11:fingers, 21:toes
  filter(of_fxcircum!=7) %>% # pathological fracture
  group_by(shareid) %>% 
  top_n(-1,of_FxDate) %>% 
  top_n(-1,of_FxSite) %>% 
  top_n(-1,of_fxside) %>%
  top_n(-1,of_fxcircum) %>%
  mutate(fx=1)
cohortfx = merge(cohort,fx,by.x="shareid",by.y="shareid",all.x=T)
cohortfx$fx[is.na(cohortfx$fx)] = 0
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
# make tables
# FX by BB
myvars = names(cohortfx)[c(6:15,18,20,21,25,27,28,35)]
catvars = names(cohortfx)[c(6,11,13:15,18,20,21,25,27,28,35)]
tab_overall = CreateTableOne(vars = myvars, data = cohortfx, factorVars = catvars)
tab_overallMat <- print(tab_overall, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tab1 <- CreateTableOne(vars = myvars, data = cohortfx, factorVars = catvars, strata = "BB")
tab1_Mat <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(cbind(tab_overallMat,tab1_Mat), file = "FxCohort.csv")
# BMD by BB
myvars = names(cohortbmd)[c(6:15,18,20,21,25,27,28,38:40,42:45)]
catvars = names(cohortbmd)[c(6,11,13:15,18,20,21,25,27,28)]
tab_overall = CreateTableOne(vars = myvars, data = cohortbmd, factorVars = catvars)
tab_overallMat <- print(tab_overall, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tab1 <- CreateTableOne(vars = myvars, data = cohortbmd, factorVars = catvars, strata = "BB")
tab1_Mat <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(cbind(tab_overallMat,tab1_Mat), file = "BMDCohort.csv")
# BMD models
pheno = c("f8cbnbmd","f8cbtobmd","f8cbtrbmd","s8cbl2bd","s8cbl3bd","s8cbl4bd","s8cbl24bd")
colnum = which(names(cohortbmd) %in% pheno)
resultsf = array(0,dim=c(length(pheno),4,2))
cohortbmdf = cohortbmd[cohortbmd$SEX==2,]
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmdf[,colnum[i]]~BB,data=cohortbmdf)
  resultsf[i,1,1] = summary(mod1)$coef[2,1]
  resultsf[i,2:3,1] = confint(mod1)[2,]
  resultsf[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmdf[,colnum[i]]~BB+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+EST8+menov+priorcvd,data=cohortbmdf)
  resultsf[i,1,2] = summary(mod2)$coef[2,1]
  resultsf[i,2:3,2] = confint(mod2)[2,]
  resultsf[i,4,2] = summary(mod2)$coef[2,4]
}
resultsm = array(0,dim=c(length(pheno),4,2))
cohortbmdm = cohortbmd[cohortbmd$SEX==1,]
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmdm[,colnum[i]]~BB,data=cohortbmdm)
  resultsm[i,1,1] = summary(mod1)$coef[2,1]
  resultsm[i,2:3,1] = confint(mod1)[2,]
  resultsm[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmdm[,colnum[i]]~BB+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmdm)
  resultsm[i,1,2] = summary(mod2)$coef[2,1]
  resultsm[i,2:3,2] = confint(mod2)[2,]
  resultsm[i,4,2] = summary(mod2)$coef[2,4]
}
resultsa = array(0,dim=c(length(pheno),4,2))
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmd[,colnum[i]]~BB,data=cohortbmd)
  resultsa[i,1,1] = summary(mod1)$coef[2,1]
  resultsa[i,2:3,1] = confint(mod1)[2,]
  resultsa[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmd[,colnum[i]]~BB+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmd)
  resultsa[i,1,2] = summary(mod2)$coef[2,1]
  resultsa[i,2:3,2] = confint(mod2)[2,]
  resultsa[i,4,2] = summary(mod2)$coef[2,4]
}
df = data.frame(pheno,cbind(resultsf[,,1],resultsf[,,2],resultsm[,,1],resultsm[,,2],resultsa[,,1],resultsa[,,2]))
lab = c("Est","Lower","Upper","p-val")
names = c(paste("Female Raw",lab),paste("Female Adj",lab))
names = c(names,paste("Male Raw",lab),paste("Male Adj",lab))
names = c(names,paste("All Raw",lab),paste("All Adj",lab))
names(df) = c("BMD Location",names)
write.table(df,file="BMDResults.txt",row.names=F,col.names=T,quote=F,sep="\t")
# FX models
mod1 = glm(fx~BB,data=cohortfx,family="binomial")
summary(mod1)
mod2 = glm(fx~BB+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortfx,family="binomial")
summary(mod2)
# BMD models for B1
pheno = c("f8cbnbmd","f8cbtobmd","f8cbtrbmd","s8cbl2bd","s8cbl3bd","s8cbl4bd","s8cbl24bd")
colnum = which(names(cohortbmd) %in% pheno)
resultsf = array(0,dim=c(length(pheno),4,2))
cohortbmdf = cohortbmd[cohortbmd$SEX==2,]
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmdf[,colnum[i]]~B1,data=cohortbmdf)
  resultsf[i,1,1] = summary(mod1)$coef[2,1]
  resultsf[i,2:3,1] = confint(mod1)[2,]
  resultsf[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmdf[,colnum[i]]~B1+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+EST8+menov+priorcvd,data=cohortbmdf)
  resultsf[i,1,2] = summary(mod2)$coef[2,1]
  resultsf[i,2:3,2] = confint(mod2)[2,]
  resultsf[i,4,2] = summary(mod2)$coef[2,4]
}
resultsm = array(0,dim=c(length(pheno),4,2))
cohortbmdm = cohortbmd[cohortbmd$SEX==1,]
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmdm[,colnum[i]]~B1,data=cohortbmdm)
  resultsm[i,1,1] = summary(mod1)$coef[2,1]
  resultsm[i,2:3,1] = confint(mod1)[2,]
  resultsm[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmdm[,colnum[i]]~B1+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmdm)
  resultsm[i,1,2] = summary(mod2)$coef[2,1]
  resultsm[i,2:3,2] = confint(mod2)[2,]
  resultsm[i,4,2] = summary(mod2)$coef[2,4]
}
resultsa = array(0,dim=c(length(pheno),4,2))
for (i in 1:length(pheno)) {
  mod1 = lm(cohortbmd[,colnum[i]]~B1,data=cohortbmd)
  resultsa[i,1,1] = summary(mod1)$coef[2,1]
  resultsa[i,2:3,1] = confint(mod1)[2,]
  resultsa[i,4,1] = summary(mod1)$coef[2,4]
  mod2 = lm(cohortbmd[,colnum[i]]~B1+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortbmd)
  resultsa[i,1,2] = summary(mod2)$coef[2,1]
  resultsa[i,2:3,2] = confint(mod2)[2,]
  resultsa[i,4,2] = summary(mod2)$coef[2,4]
}
df = data.frame(pheno,cbind(resultsf[,,1],resultsf[,,2],resultsm[,,1],resultsm[,,2],resultsa[,,1],resultsa[,,2]))
lab = c("Est","Lower","Upper","p-val")
names = c(paste("Female Raw",lab),paste("Female Adj",lab))
names = c(names,paste("Male Raw",lab),paste("Male Adj",lab))
names = c(names,paste("All Raw",lab),paste("All Adj",lab))
names(df) = c("BMD Location",names)
write.table(df,file="BMDResultsB1.txt",row.names=F,col.names=T,quote=F,sep="\t")
# FX models
resultsfx = array(0,dim=c(4,4))
mod1 = glm(fx~BB,data=cohortfx,family="binomial")
coef = coef(summary(mod1))
resultsfx[1,4] = coef[2,4]
resultsfx[1,1] = exp(coef)[2,1]
resultsfx[1,2:3] = exp(confint(mod1))[2,]
mod2 = glm(fx~BB+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortfx,family="binomial")
coef = coef(summary(mod2))
resultsfx[2,4] = coef[2,4]
resultsfx[2,1] = exp(coef)[2,1]
resultsfx[2,2:3] = exp(confint(mod2))[2,]
# FX models B1
mod1 = glm(fx~B1,data=cohortfx,family="binomial")
coef = coef(summary(mod1))
resultsfx[3,4] = coef[2,4]
resultsfx[3,1] = exp(coef)[2,1]
resultsfx[3,2:3] = exp(confint(mod1))[2,]
mod2 = glm(fx~B1+AGE8+HGT8+WGT8+CURRSMK8+DMRX8+HRX8+LIPRX8+priorcvd,data=cohortfx,family="binomial")
coef = coef(summary(mod2))
resultsfx[4,4] = coef[2,4]
resultsfx[4,1] = exp(coef)[2,1]
resultsfx[4,2:3] = exp(confint(mod2))[2,]
label = c("BB raw","BB adj","B1 raw","B1 adj")
dffx = data.frame(label,resultsfx)
names(dffx) = c("Model","Est","Lower","Upper","p-value")
write.table(dffx,file="ResultsFX.txt",row.names=F,col.names=T,quote=F,sep="\t")