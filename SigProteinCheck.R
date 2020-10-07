library(limma)
library(xlsx)
library(ggplot2)
library(readxl)
library(edgeR)
library(ggplot2)
library(broom)
library(dplyr)
library(STRINGdb)
library(igraph)
library(gtools)
library(matrixStats)
library(RCy3)

windowsFonts("Source Sans Pro" = windowsFont("Source Sans Pro"))

#Reading in and formatting data
proteinData = read_xlsx("Twins_SomaLogic_Array.xlsx")
proteinData = proteinData[, !names(proteinData) %in% c("7071", "7072")]
patientData = read_xlsx("Twins_Sample_Info.xlsx")
patientData = patientData[which(patientData$Age != 5),] ##Remove 5 year olds, too young
protSymb = read_xlsx("ProtSymbols.xlsx")
extPPI = read_xlsx("ExtendedPPI2.xlsx")

rownames(proteinData)<- proteinData$ID
proteinData$ID <- NULL
proteinDataNorm <- as.data.frame(t(proteinData))

##Normalization and transformation
proteinDataNorm = log2(proteinDataNorm)
proteinDataNorm = as.data.frame(normalizeQuantiles(proteinDataNorm))

##Saving Asthmatic type of patients
AsthmaType <- patientData$Asthma

##Saving Asthma group data
group <- patientData$"Asthma group"

##Saving age data
age <- patientData$Age
age <- ifelse(age <= 18, "11-18", ifelse(age <= 50, "19-50", "51-77"))

##Saving sex data
sex <- patientData$Sex

##Saving FEV Ratio data
FEVRatio <- patientData$Ratio
FEVRatio[FEVRatio == "."] <- NA
FEVRatio <- ifelse(FEVRatio < 0.7, "<0.7", ">=0.7")
# FEVRatio <- ifelse(FEVRatio[,1] < 0.7, "<0.7", ifelse(FEVRatio[,1] < 0.85, ifelse(FEVRatio[,2] < 15, "<0.85", ">0.85"), ">0.7"))

##Saving FeNo data
FeNo <- as.numeric(patientData$FeNo)

##Adding grouping and covariates to protein data
proteinDataNorm <- cbind(proteinDataNorm, AsthmaType, sex, age, FEVRatio, FeNo, group)
# proteinDataNorm <- data.frame(t(proteinData))
# proteinDataNorm <- proteinDataNorm[complete.cases(proteinDataNorm), ]
# proteinDataNorm <- t(proteinDataNorm)

##Testing significant proteins for asthma vs. no asthma
sigProteinAsthma = data.frame()

for (i in 1:(ncol(proteinDataNorm[,1:1129])))
{
  results <- tidy(t.test(proteinDataNorm[, i] ~ AsthmaType, data = proteinDataNorm))
  sigProteinAsthma <- rbind(sigProteinAsthma, results)
}

sigProteinAsthma <- cbind(colnames(proteinDataNorm[,1:1129]), sigProteinAsthma)
sigProteinAsthma$adjPValue <- p.adjust(sigProteinAsthma$p.value, "BH") ## P-value correction

##Testing significant proteins for asthma vs. non-asthma patients given significant covariates
sigProteinCov = data.frame()

for (i in 1:(ncol(proteinDataNorm[, 1:1129])))
{
  results <-tidy(glm(proteinDataNorm[, i] ~ AsthmaType + FEVRatio + FeNo, data = proteinDataNorm))$p.value[2]
  sigProteinCov = rbind(sigProteinCov, results)
}

sigProteinCov <- cbind(colnames(proteinDataNorm[,1:1129]), sigProteinCov)
colnames(sigProteinCov) <- c("Name", "p.value")

sigProteinCov$adjPVal <- p.adjust(sigProteinCov$p.value, "BH") ##P-value correction

##Calculating fold change
asthmaProt <- proteinDataNorm[AsthmaType == "yes", 1:1129]
nonAsthmaProt <- proteinDataNorm[AsthmaType == "no", 1:1129]

FCAsthma <- as.data.frame(foldchange(colMeans(asthmaProt), colMeans(nonAsthmaProt)))
sigProteinCov <- cbind(sigProteinCov, FCAsthma)
colnames(sigProteinCov)[4] <- "FC"

ggplot(sigProteinCov, aes(x = FC, y = -log10(p.value)) ) +
  geom_point()


##Graphing p-value of significant proteins
sigProteinCov <- sigProteinCov[which(sigProteinCov$adjPVal< 0.2),]
sigProteinCov$symb <- c("IGFBP6", "MMP8", "PRSS3P2", "CCL11")

ggplot(sigProteinCov[,], aes(x = reorder(sigProteinCov[,]$symb, sigProteinCov[,]$"p.value"), y = -log10(sigProteinCov[,]$"p.value"))) +
  geom_boxplot() +
  theme_minimal(base_size = 15, base_family = "Source Sans Pro") +
  labs(x = "Protein Symbol", y = "-log10(p-value)")

ggsave(sigProteinGraph, filename = "sigProteinGraph.tiff", device = "tiff", dpi = 200)
##Testing differences between discordant twins
sigProteinTwin = vector()

Dis = proteinDataNorm[which(group == "dis"),]
# NonAsthmaticDis = proteinDataNorm[which(AsthmaType == "no" & group == "dis"),]

AsthmaticCon = proteinDataNorm[which(group == "con"),]
NonAsthmaticNone = proteinDataNorm[which(group == "none"),]


for (i in 1:1129)
{
  results <-tidy(glm(Dis[, i] ~ AsthmaType + FEVRatio + FeNo, data = Dis))$p.value[2]
  sigProteinTwin[i] = results
}

sigProteinTwin = data.frame(sigProteinTwin)
sigProteinTwin <- cbind(colnames(proteinDataNorm[,1:1129]), sigProteinTwin)
colnames(sigProteinTwin) <- c("Name", "p.value")
sigProteinTwin$adjPVal <- p.adjust(sigProteinTwin$p.value, "BH") ##P-value correction

##Capitalizing protein symbols in protSymb
for (i in 1:ncol(protSymb))
{
  protSymb$symbolCap <- lapply(protSymb$symbol, function(v) {
    if (is.character(v))
    {
      return (toupper(v))
    }
    else
    {
      return (v)
    }
  })
}

##checking if proteins in extended PPI are in protein list
notIn = c()
count = 1
for (i in 1:nrow(extPPI))
{
  for (j in 1:2)
  {
    if(!extPPI[i, j] %in% protSymb$symbolCap)
    {
      notIn[count] <- extPPI[i, j]
      count = count + 1
    }
  }
}

notIn = notIn[!duplicated(notIn)]

# notIn = c("DEFA4", "IGF2", "CXCR4", "CXCR3", "IGFALS", "REG3A", "SERPINB3", "REG3G", "STAT3", "TCN1", "STAT6", "PPP3CC")

notIn = c("IGFALS", "CCR3", "CCR5", "CCL4", "DEFA4", "DEFA3", "CAMP", "IGF2", "VCAN", "C3AR1", "TNR", "CCR2", "ARRB1", "CXCR3", "CXCR4", "GNAI1","GNAI2" ,"GNAQ" ,"CXCL9", "GNAI3")

##Save proteins that were in  both extended PPI and in protein list into actual PPI
PPI = extPPI[which(!extPPI$"#node1" %in% notIn & !extPPI$node2 %in% notIn),]
PPI = PPI[which(PPI$combined_score >= 0.9),]
write.xlsx(PPI, "PPI.xlsx", sheetName = "realPPI")

##Save list of proteins in network into data frame
networkProt <- c(PPI$"#node1", PPI$node2)
networkProt <- as.data.frame(networkProt[!duplicated(networkProt)])
names(networkProt)[1] <- "symbol"

##Getting full name of proteins
networkProt$fullName = NA
for (i in 1:nrow(networkProt))
{
  if (networkProt$symbol[i] %in% protSymb$symbolCap)
  {
    networkProt$fullName[i] <- protSymb[which(networkProt$symbol[i] == protSymb$symbolCap),]$"...1"
  }
}

# networkProt[which(networkProt$symbol == "CCR3"),]$fullName <- "C-C motif chemokine 3"
# networkProt[which(networkProt$symbol == "CCR5"),]$fullName <- "C-C motif chemokine 5"
# networkProt[which(networkProt$symbol == "PRSS3P2"),]$fullName <- "Trypsin-2"
# networkProt[which(networkProt$symbol == "CCR2"),]$fullName <- "C-C motif chemokine 2"
# networkProt[which(networkProt$symbol == "INS"),]$fullName <- "Insulin"

networkProt[which(networkProt$symbol == "MMP1"),]$fullName <- "Interstitial collagenase"
networkProt[which(networkProt$symbol == "PRSS3P2"),]$fullName <- "Trypsin-2"
networkProt[which(networkProt$symbol == "ADRBK1"),]$fullName <- "beta-adrenergic receptor kinase 1"



##changing from factor to character
networkProt$symbol <- as.character(networkProt$symbol)

##Calcuating fold change between means of selected proteins in asthma vs. non-asthma
asthmaProt <- proteinDataNorm[AsthmaType == "yes", names(proteinDataNorm) %in% networkProt$fullName]
nonAsthmaProt <- proteinDataNorm[AsthmaType == "no", names(proteinDataNorm) %in% networkProt$fullName]

FCAsthma <- as.data.frame(foldchange(colMeans(asthmaProt), colMeans(nonAsthmaProt)))
##adding symbols to FCAsthma
FCAsthma$symbol = NA
for (i in 1:nrow(FCAsthma))
{
  FCAsthma$symbol[i] <- networkProt[which(rownames(FCAsthma)[i] == networkProt$fullName),]$symbol
}

##changing column name of FCAsthma
names(FCAsthma)[1] <- "foldChange"
##Calculating fold change between means of selected proteins in discordant asthma vs. discordant non-asthma
# AsthmaticDis = AsthmaticDis[, names(AsthmaticDis) %in% networkProt$fullName]
# NonAsthmaticDis = NonAsthmaticDis[, names(NonAsthmaticDis) %in% networkProt$fullName]
#
# FCDis <- foldchange(colMeans(AsthmaticDis), colMeans(NonAsthmaticDis))
#
# ##Calculating fold change between means of selected proteins in concordant asthma vs. no asthma
# AsthmaticCon = AsthmaticCon[, names(AsthmaticCon) %in% networkProt$fullName]
# NonAsthmaticNone = NonAsthmaticNone[, names(NonAsthmaticNone) %in% networkProt$fullName]
#
# FCNoneCon <- foldchange(colMeans(AsthmaticCon), colMeans(NonAsthmaticNone))

##saving protein interactions in list
proteins = c();
count = 1
for(i in 1:nrow(PPI))
{
  for (j in 1:2)
  {
    proteins[count] = PPI[i, j]
    count = count + 1
  }
}

##Creating PPI with igraph
asthmaFCNetwork = graph(edges = as.character(proteins))
asthmaIntNetwork = graph(edges = as.character(proteins))
nonAsthmaIntNetwork = graph(edges = as.character(proteins))
plot(asthmaFCNetwork)

##Adding fold changes as an attribute to vertices
V(asthmaFCNetwork)$FC = NA
for (i in 1:length(V(asthmaFCNetwork)))
{
  V(asthmaFCNetwork)$FC[i] <- FCAsthma[which(V(asthmaFCNetwork)$name[i] == FCAsthma$symbol),]$foldChange
}
E(asthmaFCNetwork)$score = PPI$"combined_score"
V(asthmaFCNetwork)$degree = degree(asthmaFCNetwork)

##Adding intensity expressions of proteins into asthmatic protein network
asthmaProtMeans = as.data.frame(colMeans(asthmaProt))
asthmaProtMeans = cbind(asthmaProtMeans, FCAsthma$symbol)

V(asthmaIntNetwork)$intensity = NA
for (i in 1:length(V(asthmaIntNetwork)))
{
  V(asthmaIntNetwork)$intensity[i] <- asthmaProtMeans[which(V(asthmaIntNetwork)$name[i] == asthmaProtMeans$"FCAsthma$symbol"),]$"colMeans(asthmaProt)"
}

V(asthmaIntNetwork)$degree <- degree(asthmaIntNetwork)
E(asthmaIntNetwork)$score = PPI$"combined_score"

##Adding intensity expressions of proteins into non-asthmatic protein network
nonAsthmaProtMeans = as.data.frame(colMeans(nonAsthmaProt))
nonAsthmaProtMeans = cbind(nonAsthmaProtMeans, FCAsthma$symbol)

V(nonAsthmaIntNetwork)$intensity = NA
for (i in 1:length(V(nonAsthmaIntNetwork)))
{
  V(nonAsthmaIntNetwork)$intensity[i] <- nonAsthmaProtMeans[which(V(nonAsthmaIntNetwork)$name[i] == nonAsthmaProtMeans$"FCAsthma$symbol"),]$"colMeans(nonAsthmaProt)"
}

V(nonAsthmaIntNetwork)$degree <- degree(nonAsthmaIntNetwork)
E(nonAsthmaIntNetwork)$score = PPI$"combined_score"
# rbPalPos <- colorRampPalette(c("red", "blue"))
# rbPalNeg <- colorRampPalette(c("white", "grey"))
#
#
# V(g1)$color = NA
#
# V(g1)[which(V(g1)$FC > 0)]$color <- rbPalPos(10)[as.numeric(cut(V(g1)$FC[which(V(g1)$FC > 0)], breaks = 5))]
# V(g1)[which(V(g1)$FC < 0)]$color <- rbPalNeg(10)[as.numeric(cut(V(g1)$FC[which(V(g1)$FC < 0)], breaks = 5))]
#
# layout = layout_with_fr(g1)
plot(asthmaFCNetwork, layout = layout.fruchterman.reingold(asthmaFCNetwork))

cytoscapePing()
createNetworkFromIgraph(asthmaFCNetwork,"PPI")
createNetworkFromIgraph(asthmaIntNetwork,"AsthmaPPI")
createNetworkFromIgraph(nonAsthmaIntNetwork,"NonAsthmaPPI")

##Checking mean expressions between twins
proteinDataNorm.sig <- proteinDataNorm[, names(proteinDataNorm) %in% networkProt$fullName]
proteinDataNorm.sig$twinID <- patientData$"Twin Id"

proteinDataNorm.sigMelt <- melt(proteinDataNorm.sig)



