# Created: Jun 26, 2018
# 
# Sandra Andorf, Nadeau Lab, Stanford University
###############################################################################

##This file includes tasks to get a feeling for:
## - Read in a data table
## - clean up the data
## - do some checks if the data has obvious inconsistencies


##Read in the example data.
toyData <- read.table("ToyData.csv", header=TRUE, sep=",", na.strings=c("NA","na"))
?read.table

##read.table() is the general function to read in different file types.
##There is a function to read in CSV files. It has the parameters defaulted to what is necessary for such files.
toyData1 <- read.csv("ToyData.csv")
##toyData and toyData1 are the same



##Task 1:
##How many rows are in the read in file?
##Always make sure that the row number is correct after read in.
row = nrow(toyData1)



##Task 2:
##List all the column names
##Always check if all columns are there, and especially if not some empty ones are added at the end which can happen in csv files.
colNames = colnames(toyData1)

##Task 3:
##How many case and control samples are in the file?
##Many different ways to get to the answer.
cases = toyData[which(toyData$Group == "case"),]
controls = toyData[which(toyData$Group == "control"),]


##Task 4:
##What is the gender distribution (how many males, females?)
fem = 0
male = 0

for (i in toyData$sex)
{
  if (i == "male" || i == "m")
  {
    male = male + 1
  }
  else
  {
    fem = fem + 1
  }
}

##Task 5:
##Sample.ID is a mix of the run (A or B) and the participant ID (the number after "A/B-")
##We need a column with participant IDs and one with runs. Create these.
strsplit(as.character(toyData$Sample.ID), split="-", fixed=TRUE)
##One option is to do these in a loop:
toyData <- cbind(toyData, run=NA, PID=NA)
for(i in 1:nrow(toyData)){
	Sample.ID.tmp <- strsplit( as.character(toyData[i,"Sample.ID"]), split="-" )
	##Extract the first position of the vector and make it a new column called run
	toyData[i, "run"] <- Sample.ID.tmp[[1]][1]
	##Do the same for subject
	toyData[i, "PID"] <- Sample.ID.tmp[[1]][2]
}
##I left the solution in since we didn't cover loops yet. Please try to understand this.


##Task 6:
##Check the Group asssignment for each runs of each participant.
##What do you see?


##Task 7:
##Do the same table for "run".
toyData$run
##Extract the participant IDs which have not exactly run A.
runB = toyData[which(toyData$run == "B"), ]
runB$PID
##Extract the participant IDs which have not exactly 1 run B.
runA = toyData[which(toyData$run == "A"), ]
runA$PID
##An additional check that I would do is to check if everyone has exactly 2 runs. Check that.
bool = TRUE
for (i in runB$PID)
{
  if (!(i %in% runA$PID))
  {
    bool = FALSE
  }
}

##Task 8:
##Exclude "The worst" study from toyData
toyData = toyData[!(toyData$Study.Name == "The worst"), ]

##Task 9:
##Take out all rows that are Drugs = NA and Infection_viral is no
toyData = toyData[!(is.na(toyData$Drugs) & toyData$Infection_viral == "no"), ]
##Task 10:
##How man case and control subjects are listed in the table?
cases = toyData[which(toyData$Group == "case"),]
controls = toyData[which(toyData$Group == "control"),]
##Is everthing correct?


##Task 11:
###Calculate the average age. (Do this across all rows, even though it represents several runs per person and
##normally we would build the average age over unique subjects and not the runs. But it is just to
##play a bit with data hete).

avg = mean(as.numeric(as.character(toyData$Age)), na.rm = TRUE)


##Task 12:
##Does the age vary for participants between runs?
avgA = mean(as.numeric(as.character(runA$Age)), na.rm = TRUE)
avgB = mean(as.numeric(as.character(runB$Age)), na.rm = TRUE)

