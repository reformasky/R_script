setwd("D:/thesis")

source("./R_script/util.R")
fileNamesList = c("pam50", "miller_P53DLDA");
sampleTypesList = list(c("Basal","Her2","LumA","LumB"),c("X0","X1"));
numOfCLustersList = c(4,2);
normalizationList = c("normalizationZScore", "normalizationLinear")

for(i in 1 : length(fileNamesList)) {
	sampleTypes = unlist(sampleTypesList[i]);
	sourceFileName = fileNamesList[i];

	sData = read.csv(sourceFileName, sep = "\t", header=TRUE);
	hInfo = colnames(sData);
	hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));
	sourceData = sData[,hInfo > 0 & hInfo <= length(sampleTypes)];
	sourceData = as.matrix(sourceData);
	colnames(sourceData) = hInfo[hInfo > 0 & hInfo <= length(sampleTypes)];
	
	for(j in 1 : length(normalizationList)) {
		titleName = paste(c(fileNamesList[i], normalizationList[j]), collapse = "_");
		evaluateDistanceAndClustering(sourceData, scalingMethod = normalizationList[j],numOfClusters = numOfCLustersList[i],titleName = titleName)
	}
		

}





