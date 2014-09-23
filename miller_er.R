setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./miller_er", sep = "\t", header=TRUE);
sampleTypes = c("ER");
hInfo = colnames(sData);
hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));


sourceData = sData[,hInfo > 0 & hInfo < 3];
sourceData = as.matrix(sourceData);
benchmark = hInfo[hInfo > 0 & hInfo < 3];
colnames(sourceData) = benchmark
basePath = "./plots/miller_er"
states  = 3  : 5
numOfClusters = 2
normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
for(i in 1 : length(discretization)) {
	titleName = discretization[i];
	bLine = baseLine(sourceData, numOfClusters = numOfClusters, normalization = match.fun(normalization[i]));
	similarity = pairWiseComparision(sourceData, numOfClusters = numOfClusters, benchMark = c(),
		states = states,
		normalization = normalization[i], discretization = match.fun(discretization[i]), dendro = F)
	plotEvaluations(similarity, states = states, baseLine = bLine,
	 	titleName = titleName, savePath = basePath, lx = 4.2, ly = 1)
}