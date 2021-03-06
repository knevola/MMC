vcftotidy <- function(files){
  for (i in files){
    print(i)
    gene = read.vcfR(paste(i, ".recode.vcf",sep = ""))
    # ID numbers are not unique so set ID to arbitrary numbers
    gene@fix[,3] <-1:length(gene@fix[,3])
    
    # Convert VcrF to tidy data
    tidy<- vcfR2tidy(gene, format_fields = c("GT", "DS", "GP"))
    
    # Calculating the minor allele frequency
    maf <- as.data.frame(maf(gene))
    tidy$meta
    tidy$fix
    tidy$gt
    tidy$fix$MAF <- maf$Frequency
    
    
    save(tidy, file = paste(i,"_vcfr_tidy.RData",sep=""))
    
    # Filter for MAF > 0.05 (5 SNPs)
    tidy$fix <- tidy$fix[tidy$fix$MAF > 0.05,]
    tidy$gt <- tidy$gt[tidy$gt$POS %in% tidy$fix$POS,]
    
    save(tidy, file = paste(i,"_vcfr_tidy_filtered.RData",sep=""))
    
    tidy$meta
    tidy$gt$Genotype <-NA
    tidy$gt$Genotype[tidy$gt$gt_GT == "0|0"] <- "normal"
    tidy$gt$Genotype[tidy$gt$gt_GT == "1|0" |tidy$gt$gt_GT == "0|1"] <- "heterozygous"
    tidy$gt$Genotype[tidy$gt$gt_GT == "1|1"] <- "homozygous alternative"
    tidy$gt$Genotype <- as.factor(tidy$gt$Genotype)
    table(tidy$gt$Genotype)
    save(tidy, file = paste(i,"_vcfr_tidy_filtered_genotype.RData",sep=""))
  }
}
