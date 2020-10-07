library(limma)
library(ggplot2)
library(readxl)
library(edgeR)
library(ggplot2)
proteinData = read_xlsx("Proteomics_annotated.xlsx")
# proteinData = proteinData[, !names(proteinData) %in% c("7071", "7072")]
patientData = read.csv("Samples_Phenotype.csv")
# patientData = patientData[which(patientData$Age != 5),] ##Remove 5 year olds, too young

##Formatting data
patientData = patientData[which(!colnames(patientData) %in% c("twin_id", "asthma_ever", "asthma_current", "weight", "bmi", "zygosity",
                                                              "Tech", "Twin_pair"))]
rownames(patientData) <- patientData$sampleid
proteinData = proteinData[which(!proteinData$ID %in% c("Phospholipase A2, membrane associated", "Fatty acid-binding protein, epidermal",
                                                      "Granulocyte-macrophage colony-stimulating factor",
                                                      "Vascular endothelial growth factor A, isoform 121")),]

for (i in 1:nrow(proteinData))
{
  proteinData$Symbol <- lapply(proteinData$Symbol, function(v) {
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

symbols <- proteinData$Symbol

age = patientData$age
age <- ifelse(age <= 18, "0_18", ifelse(age <= 50, "19_50", "51_77"))
age <- as.factor(age)


patientData$Asthma_Status <- as.factor(patientData$Asthma_Status)
patientData$sex <- as.factor(patientData$sex)
patientData$Gina_score <- as.factor(patientData$Gina_score)
patientData$Asthma_Concordance <- as.factor(patientData$Asthma_Concordance)

sex = patientData$sex
gina = patientData$Gina_score

proteinData = proteinData[,colnames(proteinData) %in% patientData$sampleid]
rownames(proteinData) <- symbols

##Transformation + Normalization
proteinNorm = log2(proteinData)
# proteinNorm = normalizeQuantiles(as.matrix(proteinNorm))

Asthma <- proteinNorm[, colnames(proteinNorm) %in% patientData$sampleid[patientData$Asthma_Status == "Asthmatic"]]
NonAsthma <- proteinNorm[, colnames(proteinNorm) %in% patientData$sampleid[patientData$Asthma_Status == "NonAsthmatic"]]

ks.out <- data.frame()
ks.output <- data.frame()

for (i in 1:nrow(Asthma[1:826, ])) {
  ks <- ks.test(Asthma[i,], Nonasthma[i,], alternative = "two.sided")
  ks.out <- data.frame(Gene = rownames(Asthma)[i], statistic = ks$statistic, Pvalue = ks$p.value)
  ks.output <- rbind(ks.output, ks.out)
}

ks.output$adjP <- p.adjust(ks.output$Pvalue, "fdr")






















##Saving asthma status
patientData$Asthma_Status = (ifelse(patientData$Asthma_Status == "1", "Asthmatic", "NonAsthmatic"))

Asthma_Status = patientData$Asthma_Status

##Saving twin type
patientData$Asthma_Concordance = (ifelse(patientData$Asthma_Concordance == "1", "Con", "Dis"))

Asthma_Con = patientData$Asthma_Concordance

design <- model.matrix(~0 + Asthma_Con + age + sex + gina)

# colnames(design) <- c("NonAsthmatic", "Asthmatic")
# cont.matrix = makeContrasts(Asthma_Status = Asthma_StatusAsthmatic - Asthma_StatusNonAsthmatic, levels = design)
cont.matrix = makeContrasts(Asthma_Con = Asthma_ConCon - Asthma_ConDis, levels = design)

fit <- lmFit(proteinNorm, design)

fit.cont <- contrasts.fit(fit, cont.matrix)

fit.cont <- eBayes(fit.cont, trend = T)

dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

volcanoplot(fit.cont)
topTable(fit.cont, coef="Asthma_Con", sort.by = "p")
