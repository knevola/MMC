# Figure GO and KEGG Enrichment 
rm(list = ls())
setwd("/home/clary@mmcf.mehealth.org/Framingham/OmicData/MMC/data")
library(VennDiagram)
library(dplyr)
library(ggplot2)

GO_19 <- read.csv("GO_miR19a.csv")
KEGG_19<-read.csv("KEGG_miR19a.csv")
GO_186<-read.csv("GO_miR186.csv")
KEGG_186 <-read.csv( "KEGG_miR186.csv")
GO_blue <- read.csv("BlueClusterGO.csv")
KEGG_blue <- read.csv("BlueClusterKEGG.csv")
GO_terms <- read.csv("GOterms.csv", stringsAsFactors = F)
GO_blue_all <-read.csv("BlueClusterunfilterGO.csv")
GO_19_all <- read.csv("GO_miR19a_unfiltered.csv")
GO_186_all <- read.csv("GO_miR186_unfiltered.csv")

GO_list <- list(miR_19a = GO_19$Term, miR_186 = GO_186$Term, `Blue Cluster` = GO_blue$Term)
venn.diagram(GO_list, filename = "GO_term_Overlap.png",fill = c("red", "red4", "steelblue"))

KEGG_list <- list(miR_19a = KEGG_19$Pathway, miR_186 = KEGG_186$Pathway, `Blue Cluster` = KEGG_blue$Pathway)
venn.diagram(KEGG_list, filename = "KEGG_term_Overlap.png",fill = c("red", "red4", "steelblue"))

y<-intersect(intersect(KEGG_19$Pathway, KEGG_186$Pathway), KEGG_blue$Pathway)
x<-intersect(intersect(GO_19$Term, GO_186$Term), GO_blue$Term)
KEGG_19a <- data.frame(KEGG_19, `-log10(p)` = -log10(KEGG_19$P.DE), Group = "miR-19a")  
KEGG_19a <- KEGG_19a %>% filter(., Pathway %in% y[1:5])
KEGG_186a <- data.frame(KEGG_186, `-log10(p)` = -log10(KEGG_186$P.DE), Group = "miR-186")  
KEGG_186a <- KEGG_186a %>% filter(., Pathway %in% y[1:5])
KEGG_bluea <- data.frame(KEGG_blue, `-log10(p)` = -log10(KEGG_blue$P.DE), Group = "Blue Cluster")  
KEGG_bluea <- KEGG_bluea %>% filter(., Pathway %in% y[1:5])

KEGG_all <- rbind(KEGG_19a, KEGG_186a, KEGG_bluea)
KEGG_all$Pathway <- factor(x = KEGG_all$Pathway, levels = KEGG_all$Pathway[order(KEGG_all$X.log10.p.[KEGG_all$Group == "Blue Cluster"],decreasing = F)])
ggplot(data = KEGG_all, aes(x = Pathway, y = X.log10.p., fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge()) + geom_abline(slope = 0,intercept = -log10(0.05)) + 
  coord_flip() + ylab(label = "-log10(p)") + scale_fill_manual(values = c("red", "red4", "steelblue")) + ggtitle(label = "KEGG Overlapping Pathways")

GO19_unfiltered <- data.frame(GO_19_all, logp = -log10(GO_19_all$P.DE), Group = "miR-19a")
GO186_unfiltered <- data.frame(GO_186_all, logp = -log10(GO_186_all$P.DE), Group = "miR-186")
GOblue_unfiltered <- data.frame(GO_blue_all, logp = -log10(GO_blue_all$P.DE), Group = "Blue Cluster")

GO_unfiltered <- rbind(GO19_unfiltered, GO186_unfiltered, GOblue_unfiltered)

Terms_int <- GO_unfiltered %>% filter(., Term %in% c("osteoblast differentiation", "osteoclast differentiation", "bone development", "cellular response to insulin stimulus",
                                                  "cellular response to parathyroid hormone stimulus", "regulation of intracellular estrogen receptor signaling pathway"))
Insulin <- GO_unfiltered %>%  filter(., X %in% GO_terms$Insulin)
Estrogen <- GO_unfiltered %>% filter(., X %in% GO_terms$Estrogen)
TH.PTH <- GO_unfiltered %>% filter(., X %in% GO_terms$TH.PTH)
Adrenergic <- GO_unfiltered %>% filter(., X %in% GO_terms$Adrenergic)

Terms_int$Term <- factor(x = Terms_int$Term, levels = Terms_int$Term[order(Terms_int$logp[Terms_int$Group == "Blue Cluster"],decreasing = F)])

ggplot(data = Terms_int, aes(x = Term, y = logp, fill = Group)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p)") +
  ggtitle(label = "Bone Related GO Terms") + scale_fill_manual(values = c("red", "red4", "steelblue"))

ggplot(data = Insulin, aes(x = Term, y = logp, fill = Group)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p)") +
  ggtitle(label = "Insulin Signaling GO Terms") + scale_fill_manual(values = c("red", "red4", "steelblue"))

ggplot(data = Estrogen, aes(x = Term, y = logp, fill = Group)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p)") +
  ggtitle(label = "Estrogen Signaling GO Terms") + scale_fill_manual(values = c("red", "red4", "steelblue"))

ggplot(data = TH.PTH, aes(x = Term, y = logp, fill = Group)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p)") +
  ggtitle(label = "TH/PTH GO Terms") + scale_fill_manual(values = c("red", "red4", "steelblue"))

ggplot(data = Adrenergic, aes(x = Term, y = logp, fill = Group)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_abline(slope = 0, intercept = -log10(0.05)) + coord_flip() + ylab(label = "-log10(p)") +
  ggtitle(label = "Adrenergic Signaling GO Terms") + scale_fill_manual(values = c("red", "red4", "steelblue"))
