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
sData = sData[1 : nGenes, hInfo > 0 & hInfo < length(sampleTypes)];
colnames(sData) = hInfo[hInfo > 0 & hInfo < length(sampleTypes)];

#make sure all the discretizations get the same partician of genes;
genesInGroups = generateGroups(numOfGenes = nGenes, numOfGroups = nGroups);

plotBase = "./plots/oneDJury/smallPValues/pam50"
savedFile= "./processedData/oneDJury/smallPValues/pam50.csv"
normalizations = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretizations = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfStates = 3 : 5;
numOfClusters = 4;

for (i in 1 : length(discretizations)) {
	titleName = discretizations[i];
	result = evaluateOneDJury(sData = sData, genesInGroups = genesInGroups, numOfClusters = numOfClusters, normalizations[i], discretizations[i], plotBase = plotBase)
	tiff( file= file.path(plotBase, paste(c(titleName,"tiff"), collapse = "."))
		, units="in", width= 3.5, height= 1.8, res=300);
	textplot(result, cmar = 1.2, rmar = 1.2, halign = "center", valign = "center");
	title(titleName)
	dev.off();
	colnames(result) = paste(titleName, colnames(result),   sep = " ")
	if(i == 1) {
		savedCSV = data.frame(result)
	}
	else {
		savedCSV = data.frame(savedCSV, result);
	}
}

write.csv(savedCSV, file = savedFile)