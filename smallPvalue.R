#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./smallPValues.xls", sep = "\t", header = TRUE)
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader,sampleTypes = c("Basal","Her2","LumA","LumB")))
colnames(sData) = hInfo;
#sort basing on p value 
sData = sData[order(sData$"0"),]

numOfGenes = 1 : 10 * 20;

basePath = "./plots/smallPValues"
states = 3 : 5
normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
for(nGenes in numOfGenes) 
	for(i in 1 : length(discretization)) {
		titleName = paste(c(nGenes, discretization[i]), collapse="_");
		titleName = paste(c(discretization[i],"Number of genes", nGenes), collapse = "_");
		sourceData = sData[1: nGenes, hInfo > 0 & hInfo < 5];
		colnames(sourceData) = hInfo[hInfo > 0 & hInfo < 5];
		benchmark = hInfo[hInfo > 0 & hInfo < 5];
		bLine = baseLine(sourceData, numOfClusters = 4, normalization = match.fun(normalization[i]));
		results = pairWiseComparision(sourceData, numOfClusters = 4,benchMark = benchmark,states = states, 
			normalization = normalization[i], discretization = match.fun(discretization[i]), dendro = 3, 
			savePath = basePath,titleName = titleName)
		plotEvaluations(results, states = states, baseLine = bLine,
		 	titleName = titleName, savePath = basePath, lx = 4.2, ly = 0.95)
}