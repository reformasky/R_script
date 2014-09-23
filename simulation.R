
setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)
basePath = "./plots/simulation"
normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfClusters = 5
sd = 3
states = 3 :5

generateData = function(numOfFeatures = 100,numOfSamples = 40, numOfClusters = 5, stdev = 0.6) {
	result = matrix(0,ncol = numOfSamples, nrow = numOfFeatures);
	centers = matrix(0, ncol = numOfClusters, nrow = numOfFeatures);
	labels = 1 : numOfSamples;
	for(i in 1 : numOfFeatures) {
		centers[i,] = sample(1 : numOfClusters)
	}

	for ( i in 1 : numOfSamples) {
		labels[i] = sample(1 : numOfClusters,size = 1);
		result[,i] = rnorm(numOfFeatures, mean = centers[,labels[i]], stdev);
	}
	colnames(result) = labels;
	result;
}

sourceData = generateData(stdev = sd);
benchMark = as.numeric(colnames(sourceData))

for(i in 1 : length(discretization)) {
	titleName = paste(c(discretization[i],"simulation", "sd", sd), collapse = "_")
	bLine = baseLine(sourceData, numOfClusters, normalization = match.fun(normalization[i]))
	similarity = pairWiseComparision(sourceData, numOfClusters = numOfClusters, states = states, 
		normalization = normalization[i], discretization = match.fun(discretization[i]), benchMark = benchMark, 
		dendro = 3, savePath = basePath,titleName = paste(c(titleName,"dendrograph"), collapse = "_"));
	plotEvaluations(similarity, states, baseLine = bLine, 
		titleName = titleName, savePath = basePath, lx = 4.2, ly = 0.92)
}

