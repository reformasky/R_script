library(phyclust, quiet = TRUE)
setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")
set.seed(1024)


#return a list of length numOfGroups, 
generateGroups = function(numOfGenes, numOfGroups) {
	genesInGroups = sample(1 : numOfGenes);
	genesPerGroup = numOfGenes / numOfGroups;
	result = list();
	for(g in 1 : numOfGroups) {
		result[g] = list(genesInGroups[((g - 1) * genesPerGroup + 1) : (g * genesPerGroup)])
	}
	result;
}



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
trueLabel = hInfo[hInfo > 0]
#make sure all the discretizations get the same partician of genes;
genesInGroups = generateGroups(numOfGenes = nGenes, numOfGroups = nGroups);

plotBase = "./plots/oneDJury/smallPValues/pam50"
savedFile= "./processedData/oneDJury/smallPValues/pam50.csv"
normalizations = c("normalizationLinear", "normalizationLinear", "normalizationZScore")
discretizations = c("discretizationQuantile","discretizationFloor", "discretizationZScore")
numOfStates = 3 : 5;
numOfClusters = 4;

for (i in 1 : length(discretizations)) {
	# titleName = discretizations[i];
	results = evaluateOneDJury(sData = sData, genesInGroups = genesInGroups, numOfClusters = numOfClusters, normalizations[i], discretizations[i])
	# discretized = results$discretized
	# noDiscretized = results$noDiscretized;

	# plotEvaluations(results, numOfStates, baseLine = 0, titleName = titleName, savePath = plotBase)

	# colnames(discretized) = paste(titleName, colnames(discretized),   sep = " ")
	if(i == 1) {
		savedCSV = data.frame( results)
	}
	else {
		savedCSV = data.frame(savedCSV,  results);
	}

}
 write.csv(savedCSV, file = savedFile)