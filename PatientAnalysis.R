library("ggplot2")
library("xlsx")
library("readxl")
library("tidyverse")
library("extrafont")
library("reshape2")
library("ggpubr")
library("lemon")
patientData = read_xlsx("Twins_Sample_Info.xlsx")
patientData = patientData[which(patientData$Age != 5),]

##Creating data set with only wanted variables
FEVRatioData <- patientData[, c("Sample ID", "Ratio", "Asthma", "Age")]
SexData <- patientData[, c("Sample ID", "Sex", "Asthma")]
AsthmaGroupData <- patientData[, c("Sample ID", "Asthma group")]
SexAsthmaData <- cbind(SexData, AsthmaGroupData)
AllergyData <- patientData[, c("Sample ID", "Allergies Survey [y/n]", "Allergies Survey", "allergy concordance", "Asthma")]
SmokingData <- patientData[, c("Sample ID", "SHS", "Asthma")]
AgeData <- patientData[, c("Sample ID", "Age", "Asthma")]
SexData$Type = paste(SexData$Sex, SexData$Asthma, sep="_")
SexAsthmaData$SexAsthmaGroup = paste(SexAsthmaData$Sex, SexAsthmaData$"Asthma group", sep = "_")

##Cleaning data
FEVRatioData <- FEVRatioData[FEVRatioData$Ratio != ".", ]

AllergyData = AllergyData[AllergyData$"Allergies Survey [y/n]" != "NA",]
AllergyData = AllergyData[AllergyData$"Allergies Survey [y/n]" != "don't know", ]
AllergyData[which(AllergyData$"Allergies Survey [y/n]" != "yes" & AllergyData$"Allergies Survey [y/n]" != "no"),]$"Allergies Survey [y/n]" <- "yes"

##Creating age groups for age data
AgeData$ageDistr = NA
AgeData$ageDistr <- ifelse(AgeData$Age <= 18, "11-18", ifelse(AgeData$Age <= 50, "19-50", "51-77"))

##Creating age groups for FEV Ratio data
FEVRatioData$ageDistr <- ifelse(FEVRatioData$Age <= 18, "11-18", ifelse(FEVRatioData$Age <= 50, "19-50", "51-77"))

##Getting frequencies of ages/age groups
AgeData.dat = data.frame(table(AgeData$Age))
AgeData.ageDistrDat = data.frame(table(AgeData$ageDistr))

##Matching asthma with SHS
SmokingData$Asthma_SHS = paste(SmokingData$Asthma, SmokingData$SHS, sep = "_")

##Getting frequency for asthma + SHS
SmokingData = data.frame(table(SmokingData$Asthma_SHS))

##Twin A and Twin B
twinA = data.frame()
twinB = data.frame()
for (i in 1:78)
{
  if (!i %% 2)
  {
    twinB <- rbind(twinB, patientData[i,])
    next()
  }
  twinA <- rbind(twinA, patientData[i,])
}

##Saving Twin FEV1/FVC data
twinA <- cbind.data.frame(twinA$Ratio, twinA$"Twin Id", twinA$Asthma, twinA$"Asthma group")
twinB <- cbind.data.frame(twinB$Ratio, twinB$"Twin Id", twinB$Asthma, twinB$"Asthma group")

twinA$twin <- rep("twinA", 39)
twinB$twin<- rep("twinB", 39)

colnames(twinA)[1:4] <- c("Ratio", "Twin_ID", "Asthma_Status", "Concordance")
colnames(twinB)[1:4] <- c("Ratio", "Twin_ID", "Asthma_Status", "Concordance")

twinA <- twinA[twinA$Ratio != ".", ]
twinB <- twinB[twinB$Ratio != ".", ]

twinFEV <- rbind(twinA, twinB)
twinFEV$Ratio <- as.numeric(as.character(twinFEV$Ratio))
twinFEV$Color <- ifelse(twinFEV$Concordance == "con", "Concord_Asthma", ifelse(twinFEV$Concordance == "dis", "Discordant_Asthma", "Concord_NonAsthma"))
twinFEV$Asthma_Status <- ifelse(twinFEV$Asthma_Status == "yes", "Asthmatic", "Non-Asthmatic")

##Saving allergy + asthma data into one column, and getting frequency of each type
AllergyData$Asthma_Allergies = paste(AllergyData$Asthma, AllergyData$"Allergies Survey [y/n]", sep = "_")
AllergyData.distr = data.frame(table(AllergyData$Asthma_Allergies))

colnames(SmokingData)[colnames(SmokingData) == "Var1"] <- "Asthma_SHS"
colnames(AllergyData.distr)[colnames(AllergyData.distr) == "Var1"] <- "Asthma_Allergies"

for (i in AgeData.dat$Var1)
{
  ageAsthma = AgeData[which(AgeData$Age == i & AgeData$Asthma == "yes"),]
  AgeData.dat$FreqAsthma[which(AgeData.dat$Var1 == i)] = length(ageAsthma$Age)
}

for (i in AgeData.ageDistrDat$Var1)
{
  ageAsthma = AgeData[which(AgeData$ageDistr == i & AgeData$Asthma == "yes"),]
  AgeData.ageDistrDat$FreqAsthma[which(AgeData.ageDistrDat$Var1 == i)] = length(ageAsthma$ageDistr)
  AgeData.ageDistrDat$FreqNonAsthma = AgeData.ageDistrDat$Freq - AgeData.ageDistrDat$FreqAsthma
}

AgeData.ageDistrDatMelt = melt(AgeData.ageDistrDat)
AgeData.ageDistrDatMelt = AgeData.ageDistrDatMelt[which(AgeData.ageDistrDatMelt$variable != "Freq"),]
AgeData.ageDistrDatMelt = AgeData.ageDistrDatMelt[rev(rownames(AgeData.ageDistrDatMelt)),]
SexAsthmaData = data.frame(table(SexAsthmaData$SexAsthmaGroup))
AsthmaGroupData = data.frame(table(AsthmaGroupData$"Asthma group"))
SexData = data.frame(table(SexData$Type))
SexData$gend = c(rep("Female", 2), rep("Male", 2))
SexData$asthma = c("no", "yes", "no", "yes")
##Changing ratio data to numeric, so graph will treat y-values as numbers instead of characters
FEVRatioData$Ratio = as.numeric(FEVRatioData$Ratio)

##Gettings fonts
windowsFonts("Source Sans Pro" = windowsFont("Source Sans Pro"))

#Graphs
# twinTypeGraph <- ggplot(AsthmaGroupData, aes(x = Var1, y = Freq)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Distribution of Twin Types", x = "Twin Type", y = "Frequency") +
#   theme_minimal(base_size = 20, base_family = "Source Sans Pro")
# sexAsthmaGraph <- ggplot(SexAsthmaData, aes(x = Var1, y = Freq)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Distribution of Gender and Twin Type", subtitle = "in Patient Data", x = "Gender_TwinType", y = "Frequency") +
#   theme_minimal(base_size = 20, base_family = "Source Sans Pro")
gendAsthmaGraph <- ggplot(SexData, aes(fill = SexData$asthma, x = SexData$gend, y = SexData$Freq)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Distribution of asthma", subtitle = "by gender", x = "Gender", y = "Frequency") +
  theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  scale_fill_manual("Asthma Type", labels = c("Non-Asthmatic", "Asthmatic"), values = c("yes" = "#00BFC4", "no" = "#F8766D"))
ageGraph <- ggplot(AgeData.ageDistrDatMelt, aes(fill = AgeData.ageDistrDatMelt$variable, x = Var1, y = value)) +
  geom_bar(position="dodge", stat = "identity") +
  labs(title = "Age Distribution", subtitle = "in Patient Data", x = "Age Group", y = "Frequency") +
  theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  scale_fill_manual("Asthma Type", labels = c("Asthmatic", "Non-Asthmatic"), values = c("FreqAsthma" = "#00BFC4", "FreqNonAsthma" = "#F8766D"))
# ggplot(SmokingData, aes(x = Asthma_SHS, y = Freq)) +
#   geom_bar(stat = "identity")
# ggplot(AllergyData.distr, aes(x = Asthma_Allergies, y = Freq)) +
#   geom_bar(stat = "identity")
FEVAgeGraph <- ggplot(FEVRatioData, aes(x = ageDistr, y = Ratio)) +
  geom_boxplot(aes(group = ageDistr)) +
  theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.85, linetype = "dashed", color = "blue") +
  geom_jitter() +
  labs(title = "Age Group versus FEV1/FVC Ratio", x = "Age Group", y = "FEV1/FVC Ratio")
twinFEVGraph <- ggplot(twinFEV, aes(x = twin, y = Ratio, group=Twin_ID )) +
  geom_boxplot(aes(group = twin), outlier.shape = NA) +
  geom_jitter(position = position_dodge(0.2), color = ifelse(twinFEV$Asthma_Status == "yes", "green", "red")) +
  geom_line(position = position_dodge(0.2), color = "gray") +
  theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  labs(title = "FEV1/FVC Ratios across pairs of twins", x = "Twin", y ="FEV1/FVC Ratio")
ggplot(twinFEV, aes(x = Asthma_Status, y = Ratio )) +
  geom_boxplot(aes(group = Asthma_Status), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.2), aes(color = twinFEV$Color)) +
  geom_line(position = position_jitterdodge(dodge.width = 0.2), aes(group = factor(twinFEV$Twin_ID)), color = "gray") +
  # theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  labs(title = "FEV1/FVC Ratios across pairs of twins", x = "Twin", y ="FEV1/FVC Ratio")
twinFEVGraph <- ggplot(twinFEV, aes(x = Asthma_Status, y = Ratio )) +
  geom_boxplot(aes(group = Asthma_Status), outlier.shape = NA) +
  geom_pointpath(aes(color = twinFEV$Color), group = twinFEV$Twin_ID, position = position_jitter(0.2), distance = NA) +
  theme_minimal(base_size = 20, base_family = "Source Sans Pro") +
  scale_fill_manual(aesthetics = "color", "Concordance", labels = c("Concordant Asthma", "Concordant Non-Asthma", "Discordant Asthma"),
                    values = c("Concord_Asthma" = "#F8766D", "Concord_NonAsthma" = "#00BA38", "Discordant_Asthma" = "#619CFF")) +
  labs(title = "FEV1/FVC Ratios Across Pairs of Twins", x = "Asthma Status", y ="FEV1/FVC Ratio")

##Saving images
ggsave(gendAsthmaGraph, filename = "gendAsthmaGraph.tiff", device = "tiff", dpi = 200)
# ggsave(twinTypeGraph, filename = "twinTypeGraph.tiff", device = "tiff")
# ggsave(sexAsthmaGraph, filename = "sexAsthmaGraph.tiff", device = "tiff")
ggsave(ageGraph, filename = "ageGraph.tiff", device = "tiff", dpi = 200)
ggsave(FEVAgeGraph, filename = "FEVAgeGraph.tiff", device = "tiff", dpi = 200)
ggsave(twinFEVGraph, filename = "TwinFEVGraph.tiff", device = "tiff", dpi = 200)
