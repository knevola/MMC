rm(list = ls())
library(tidyverse)
library(emmeans)
library(limma)
library(WGCNA)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data")
sample_data <- read.delim("/home/clary@mmcf.mehealth.org/Framingham/OmicData/data/phe000002.v7_release_manifest.txt", comment.char="#")
data_adj <- read_delim("FinalFile_Gene_OFF_2446_Adjusted_c1.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pheno <- read.csv('PhenoData_5_28.csv', stringsAsFactors = T)

# Pheno Tscore and grouping
# pheno$Tscore <- ((-0.023 + 0.939 * pheno$f8cbnbmd  - 0.019)/1.087 - 0.858) / 0.120
# 
# 
# pheno$BMDgroup[pheno$Tscore >= -1]<- 'normal'
# pheno$BMDgroup[pheno$Tscore < -1]<- 'low'
# 
# 
# pheno$BMDgroup<- as.factor(pheno$BMDgroup)
# BMDgroup<- pheno$BMDgroup
#pheno$bbyes <- as.factor(pheno$bbyes)
#bbyes  <- pheno$bbyes
#pheno$b1yes <- as.factor(pheno$b1yes)
#b1yes <- pheno$b1yes

# merge sample data with pheno data
i<-intersect(sample_data$Subject_ID, pheno$ShareID)

pheno_samp <- merge(pheno, sample_data[,c(1,2)], by.x = "shareid", by.y = "Subject_ID")

# Merge Pheno_samp with data_adj
data_adj_t <-as.data.frame(t(data_adj))
data_adj_t <- data_adj_t[-1,]
colnames(data_adj_t)<- data_adj$transcript_cluster_id

dat<- merge(pheno_samp, data_adj_t, by.x = "Sample_ID", by.y = 0)

dataExpr <- dat[,-c(1:23)]

gag <- goodSamplesGenes(dataExpr, verbose = 0)
gag$allOK # True so dont need to remove any genes

# Variance Filtering ------------------------------------------------------
thevar <- diag(var(dataExpr))
vars<- quantile(thevar)
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/Clustering_mRNA")
write.csv(vars, "Var_mRNA.csv")

vars
qt75 <- quantile(thevar, probs = 0.75)
thebigvar <- thevar[thevar > qt75]
bigvargenes <- names(thebigvar)
row.names(dataExpr)<- dat$ShareID

dataExpr_bv <-dataExpr[,bigvargenes]
write.csv(dataExpr_bv, 'mRNA_bigvar_data.csv')

# Sample Clustering -------------------------------------------------------

sampleTree = hclust(dist(dataExpr_bv), method = "average")


# Graph Sample Clustering
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample Clustering to detect outliers", sub = "",
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Look at clustering by traits
datTraits<-dat[1:23]
datTraits1<-datTraits[-c(1:2)]
datTraits1$SEX <- as.numeric(datTraits1$SEX)
datTraits1$CURRSMK8 <- as.numeric(datTraits1$CURRSMK8)
datTraits1$DMRX8 <- as.numeric(datTraits1$DMRX8)
datTraits1$HRX8 <- as.numeric(datTraits1$HRX8)
datTraits1$LIPRX8 <- as.numeric(datTraits1$LIPRX8)
datTraits1$EST8 <- as.numeric(datTraits1$EST8)
datTraits1$BB <- as.numeric(datTraits1$BB)
datTraits1$B1 <- as.numeric(datTraits1$B1)
datTraits1$menov <- as.numeric(datTraits1$menov)
datTraits1$priorcvd <- as.numeric(datTraits1$priorcvd)

traitColors = numbers2colors(datTraits1, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits1),
                    main = "Sample dendrogram and trait heatmap")
# Sex and Estrogen clustering a bit

save(dataExpr_bv, datTraits1, file = "FHS_mRNA_1.RData")

# module detection --------------------------------------------------------
#Automatic network construction and module detection
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Call the network topology analysis function
sft = pickSoftThreshold(dataExpr_bv, powerVector = powers, verbose = 5)
#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "blue")
#This line corresponds to used an R^2 cutoff of h
abline(h=0.9, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

# Choosing a power of 4 based on mean connectivity

#Constructing the gene network and identifying modules: auto
cor <- WGCNA::cor
net = blockwiseModules(dataExpr_bv, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
#How many modules were identified? (11/10)
cor <- stats::cor
modulebkdn<- table(net$colors)
write.csv(modulebkdn, "Genes_per_module.csv")

#View heirarchical clustering
#Open a graphics window
sizeGrWindow(12,9)
#convert labels to colors for ploting
mergedColors = labels2colors(net$colors)
#Plot the dendrogram
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "mRNA-networkConstruction-auto.RData")

# Relate Modules to Traits ----------------------------------------
#Recalculate MEs with color labels
nGenes = ncol(dataExpr_bv)
nSamples = nrow(dataExpr_bv)

MEs0 = moduleEigengenes(dataExpr_bv, moduleColors)$eigengenes
row.names(MEs0)<- datTraits$shareid
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Graphical Representation
sizeGrWindow(10,6)
#Display correlations and their Pvalue
textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3,3))

#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(datTraits1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5, 
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#Names of modules (colors)
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(dataExpr_bv, MEs, use = "p"))
geneModuleMembership$color <-mergedColors
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
write.csv(geneModuleMembership, "genemodulemembership.csv" )
write.csv(MEs, "mRNA_MEs.csv")
write.csv(geneModuleMembership, "mRNA_geneModuleMembership.csv")
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
