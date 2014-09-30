library(phyclust, quiet = TRUE)
setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")
set.seed(1024)


#total number of genes for consideration
nGenes  = 500;
#numOfGroups : the length of 1DJury vector of each sample
nGroups = 20;
#all genes will be used for consideration for 1D jury score, and the distribution for 


sData = read.csv("./smallPValues.xls", sep = "\t", header = TRUE)
sampleTypes = c("Basal","Her2","LumA","LumB")
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader,sampleTypes = sampleTypes))
colnames(sData) = hInfo;
#sort basing on p value 
sData = sData[order(sData$"0"),]
#currently only consider the lowest 500 hundred genes.
sData = sData[1 : nGenes, hInfo > 0 & hInfo <= length(sampleTypes)];
colnames(sData) = hInfo[hInfo > 0 & hInfo <= length(sampleTypes)];

#make sure all the discretizations get the same partician of genes;
genesInGroups = generateGroups(numOfGenes = nGenes, numOfGroups = nGroups);

plotBase = "./plots/oneDJury/smallPValues/pam50"
savedFile= "./processedData/oneDJury/smallPValues/pam50.csv"
normalizations = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretizations = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfStates = 3 : 5;
numOfClusters = 4;

for (i in 1 : length(discretizations)) {
	i = 3
	titleName = discretizations[i];
	results = evaluateOneDJury(sData = sData, genesInGroups = genesInGroups, numOfClusters = numOfClusters, normalizations[i], discretizations[i], plotBase = plotBase)
	discretized = results$discretized
	noDiscretized = results$noDiscretized;

	plotEvaluations(results, numOfStates, baseLine = 0, titleName = titleName, savePath = plotBase)

	colnames(discretized) = paste(titleName, colnames(discretized),   sep = " ")
	if(i == 1) {
		savedCSV = data.frame(results)
	}
	else {
		savedCSV = data.frame(savedCSV, results);
	}

}
 write.csv(savedCSV, file = savedFile)