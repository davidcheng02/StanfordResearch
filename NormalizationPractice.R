library(limma)
library(ggplot2)
library(readxl)

proteinData = read_xlsx("Twins_SomaLogic_Array.xlsx")

rownames(proteinData)<- proteinData$ID
proteinData$ID <- NULL
colnames(proteinData) <- sub("X", "", colnames(proteinData))
proteinData <- data.frame(proteinData)

proteinNorm = log2(proteinData)
proteinNorm  = normalizeQuantiles(proteinNorm)

proteinNormT = data.frame(t(proteinNorm))

proteinNormT$ID <- rownames(proteinNormT)

library("reshape2")
protein.dat <- melt(proteinNormT)
protein.dat$ID <- sub("X", "", protein.dat$ID)

proteinDataGraph <- ggplot(protein.dat, aes(x = ID , y = value)) +
                            geom_boxplot()

ggsave(proteinDataGraph, filename = "proteinDataGraph.png", device = "png")
summary(proteinData)
