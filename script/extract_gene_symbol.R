# Extract Gene Symbol #
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data")

proteins <- read.csv(file = 'itraqproteins.csv')

proteins$GeneSymbol<- gsub(proteins$Description,pattern = '.*Gene Symbol = ', replacement = '')
proteins$GeneSymbol <- gsub(proteins$GeneSymbol, pattern = ').*', replacement = '')

proteins <- proteins[-1]
proteins<- unique(proteins[,1:3])

write.csv(proteins,'itraqproteins_unique.csv', row.names = F)
