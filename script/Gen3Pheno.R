#Phenotype Data: Gen 3#
# Setup R ####
rm(list = ls())
library(tidyverse)
setwd("/home/clary@mmcf.mehealth.org/Framingham")
# Read in Data ####
date <- read.delim("phs000007.v30.pht003099.v5.p11.c1.vr_dates_2014_a_0912s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)
cov <- read.delim('phs000007.v30.pht006026.v2.p11.c1.vr_wkthru_ex02_3b_0464s.HMB-IRB-MDS.txt',
                  skip=10,header=T,stringsAsFactors = F)
drug <- read.delim("phs000007.v30.pht003098.v4.p11.c1.vr_meds_2011_m_0675s.HMB-IRB-MDS.txt", 
                   skip=10,header=T,stringsAsFactors = F)
bmd <- read.delim("phs000007.v30.pht001182.v6.p11.c1.t_bmdhs_2009_3_0627s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)
bmd626 <- read.delim("phs000007.v30.pht001892.v5.p11.c1.t_bmdhs_2010_3_0626s.HMB-IRB-MDS.txt",
                     skip=10,header=T,stringsAsFactors = F)
bmd625 <- read.delim("phs000007.v30.pht002346.v1.p11.c1.t_bmdhs_2011_3_0625s.HMB-IRB-MDS.txt",
                     skip=10,header=T,stringsAsFactors = F)
est <- read.delim("phs000007.v30.pht000693.v5.p11.c1.vr_meno_ex02_3_0653s.HMB-IRB-MDS.txt",
                  skip=10,header=T,stringsAsFactors = F)
mi <- read.delim("phs000007.v30.pht003335.v6.p11.c1.vr_soe4srv_2014_a_1027s.HMB-IRB-MDS.txt",
                 skip=10,header=T,stringsAsFactors = F)
# Filter for gen 3 participants ####
gen3 <- date %>% dplyr::select(shareid,age2,att2,date2,idtype) %>% 
  filter(idtype==3) %>% # In Gen3 Cohort
  filter(att2==1) # Attended Exam 2

# Filter for cov of interest: Exam 2 ####
cov2 = cov %>% # 
  dplyr::select(shareid,SEX,AGE2,BMI2,HGT2,WGT2,CURRSMK2,CPD2,DMRX2,HRX2,LIPRX2,IDTYPE) %>% 
  filter(IDTYPE==3)

# BB Use Gen3 ####
drug2 <- drug %>% 
  dplyr::select(shareid,ther_gp1,medname,medprn,chem_gp1,chem_nm1) %>% 
  filter(ther_gp1=="BETA BLOCKING AGENTS") %>% 
  filter(!(medname %in% c("ISTALOL","IC ATENOLOL"))) %>% 
  filter(medprn==0) %>% # 5 PRN
  mutate(B1drug=ifelse(chem_gp1=="Beta blocking agents, selective",1,0)) %>% 
  mutate(BBdrug=1)  %>% 
  group_by(shareid) %>%
  summarize(BB=max(BBdrug),B1=max(B1drug))

# MI Gen3 Cohort ####
mi3 <- mi[mi$IDTYPE ==3,]

# Exam 2 Estrogen Levels ####
est2 <- est %>%
  dplyr::select(shareid,OVREM2,EST2,STOPAGE)
# BMD Data Gen3 ####
bmd2 <- rbind(bmd, bmd626)
bmd2 <- rbind(bmd2, bmd625)
bmd2 <- merge(bmd2, gen3,by.x="shareid",by.y="shareid")
bmd2 <- bmd2 %>%
  mutate(f_dpast2 = f2scdt-date2) %>% #Days femurs scan after exam 2 date
  mutate(s_dpast2 = s2scdt-date2) %>% # Days spine scan after exam 2 date
  mutate(mis_f2 = is.na(f2scdt) | f_dpast2<0) %>% # Is femur data before exam 2 or missing
  mutate(mis_s2 = is.na(s2scdt) | s_dpast2<0) %>% # Is spine data before exam 2 or missing
  filter(!mis_f2 | !mis_s2) %>% # Filter for non missing or before exam 8 femur of spine data
  dplyr::select(shareid,f2nbmd, f2tobmd, f2trbmd, s2l2bd, s2l3bd, s2l4bd, s2l24bd)

# Merge Phenotype data together ####
cohort <- merge(gen3,cov2,by.x="shareid",by.y="shareid")
cohort <- merge(cohort, drug2,by.x="shareid",by.y="shareid",all.x=T)
cohort <- merge(cohort,est2,by.x="shareid",by.y="shareid",all.x=T)
cohort <- merge(cohort, bmd2, by.x ="shareid",by.y="shareid")
cohort$BB[is.na(cohort$BB)] = 0
cohort$B1[is.na(cohort$B1)] = 0

# Create menov variable ####
cohort <- cohort %>% mutate(menov = (OVREM2==2 | STOPAGE<age2)) # 2 ovaries removed or menopause
cohort$menov[cohort$SEX==1] <- 2 # If male <- 2
cohort$EST2[cohort$SEX==1] <- 2 # If male <- 2

# Create priorcvd variable ####
cohortmi <- merge(cohort,mi3,by.x="shareid",by.y="shareid",all.x=T)
cvd <- cohortmi %>%
  mutate(cvd_event = ifelse(is.na(EVENT),0,ifelse(EVENT %in% c(1:19,30:49),1,0))) %>%
  mutate(event_days8 = DATE-date2) %>%
  mutate(priorcvd_event = ifelse(event_days8<0 & cvd_event==1,1,0)) %>%
  group_by(shareid) %>%
  summarize(priorcvd = max(priorcvd_event))
cohort <- merge(cohort, cvd,by.x="shareid",by.y="shareid",all.x=T)
cohort$priorcvd[is.na(cohort$priorcvd)] <- 0 # Set NAs to 0

write.csv(cohort, "GEN3Pheno.csv")

