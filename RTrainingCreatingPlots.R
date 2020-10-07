toyData <- read.csv("ToyData.csv")

cases = toyData[which(toyData$Group == "case"),]
controls = toyData[which(toyData$Group == "control"),]

caseSumTCell = sum(controls$T.cells, na.rm = TRUE)
controlSumTCell = sum(controls$T.cells, na.rm = TRUE)

caseSumBCell = sum(cases$B.cells, na.rm = TRUE)
controlSumBCell = sum(controls$B.cells, na.rm = TRUE)

caseSumTReg = sum(cases$Tregs, na.rm = TRUE)
controlSumTReg = sum(controls$Tregs, na.rm = TRUE)

yvalues = c(caseSumTCell, controlSumTCell, caseSumBCell, controlSumBCell, caseSumTReg, controlSumTReg)
barplot(yvalues, names.arg = c("TCell", "BCell","TReg"), col = c("darkblue", "red"), beside = T)

data <- toyData[, c("Group", "T.cells", "B.cells", "Tregs")] 

library("ggplot2")

data.trans <- melt(data)
data.trans$com <- paste(data.trans$Group, data.trans$variable, sep="_")


ggplot(data.trans, aes(x=variable, y=value)) +
  geom_violin(aes(group = data.trans$com, color = data.trans$com))