rm(list = ls())
library(GEOquery)
options(stringasFactor =F)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/MarkerInfo/phg000835.v2.FHS_SHARE_imputed_HRC1.marker-info.MULTI/phg000835.v1.FHS_SHARE_imputed_HRC1.marker-info.MULTI/snpinfo")
files<-list.files(path = "/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/MarkerInfo/phg000835.v2.FHS_SHARE_imputed_HRC1.marker-info.MULTI/phg000835.v1.FHS_SHARE_imputed_HRC1.marker-info.MULTI/snpinfo",
           pattern = "*.gz")
files <- paste("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/MarkerInfo/phg000835.v2.FHS_SHARE_imputed_HRC1.marker-info.MULTI/phg000835.v1.FHS_SHARE_imputed_HRC1.marker-info.MULTI/snpinfo/", files, sep = "")
for (i in 1:length(files)){
  print(i)
  gunzip(filename = files[i], overwrite = T, remove = F)
}

info_files <- list.files(path = "/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/MarkerInfo/phg000835.v2.FHS_SHARE_imputed_HRC1.marker-info.MULTI/phg000835.v1.FHS_SHARE_imputed_HRC1.marker-info.MULTI/snpinfo",
                         pattern = "*.info$")
chrs<- read.table(info_files[1], header = T)
for (i in 2:length(info_files)){
  print(i)
  chrs <- rbind(chrs,read.table(info_files[i], header = T))
}

chrs8 <- chrs[chrs$Rsq > 0.8,]


write.csv(chrs8, "WellImputedSNPs08.csv", quote = F, row.names = F)


setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data")
library(tidyr)
chrs8<-read.csv("WellImputedSNPs08.csv")
chrs8 <- separate(data=chrs8, col = SNP, into = c("CHR", "POS"), sep = ":")
chrs8$CHR <- as.numeric(chrs8$CHR)
chrs8$POS <- as.numeric(chrs8$POS)
chrs8pos <- chrs8[1:2]
write.table(chrs8pos,file = "WellImputedPositions8.txt", quote = F, row.names = F, sep = '\t')

vcfgzfiles<- list.files(path = "/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/Genetic Association Study/other genes/data",pattern = ".vcf.gz")

fileConn <- file("gzfile.txt")    
writeLines(vcfgzfiles, fileConn)    

unique(chrs8pos$CHR)
