library(limma)
library(ggplot2)
library(readxl)
library(edgeR)
library(ggplot2)
proteinData = read_xlsx("Twins_SomaLogic_Array.xlsx")
proteinData = proteinData[, !names(proteinData) %in% c("7071", "7072")]
patientData = read_xlsx("Twins_Sample_Info.xlsx")
patientData = patientData[which(patientData$Age != 5),] ##Remove 5 year olds, too young

rownames(proteinData)<- proteinData$ID
proteinData$ID <- NULL

##Normalize data
proteinDataNorm = log2(proteinData)
proteinDataNorm = normalizeMedianValues(as.matrix(proteinDataNorm))

##Transposing data
proteinDataNorm = as.data.frame(t(proteinDataNorm))

##Storing Asthmatic data
AsthmaType <- (ifelse(patientData$Asthma == "yes", 1, 0))
AsthmaType.fact <- factor(AsthmaType)

##Storing Age data
ages = patientData$Age

##Storing gender data
genders <- patientData$Sex

##Storing FEV da
FEVRatio <- patientData$Ratio
FEVRatio[FEVRatio == "."] <- NA
FEVRatio <- ifelse(FEVRatio < 0.7, "<0.7", ">=0.7")

##Storing SHS Data
SHS <- patientData$SHS

##Storing Smoking Data
smoke <- patientData$"ever smoker"
smoke <- ifelse(smoke == "n" | smoke == "No" | smoke == "no", "no", "yes" )

##Storing FeNo Data
FeNo <- as.numeric(patientData$FeNo)

##Storing Allergy Data
allergy <- patientData$"Allergies Survey [y/n]"
allergy <- ifelse(allergy == "don't know", "no", ifelse(allergy != "yes" & allergy != "no" & !is.na(allergy), "yes", allergy))

##Generalized linear model
fit <- glm(formula = AsthmaType ~ ages + FEVRatio + FeNo + genders + SHS + smoke + allergy, data = proteinDataNorm, family = gaussian(link = "identity"))

results <- t.test(AsthmaType ~ genders, data = proteinDataNorm)
