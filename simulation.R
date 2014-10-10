
setwd("D:/thesis")
source("./R_script/util.R")
set.seed(512)
plotPath = "./plots/simulation"
saveBase = "./processedData/discretizations/simulation"
normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfClusters = 5
sds = 3 : 5
states = 3 :5

generateData = function(numOfFeatures = 100,numOfSamples = 200, numOfClusters = 5, stdev) {
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


for(sd in sds){	
	sourceData = generateData(stdev = sd);
	plotBase = file.path(plotPath, paste(c("sd", sd), collapse = "_"))
	saveFileName =  paste(c("sd", sd), collapse = "_");
	
	plotBarGraph(sourceData, numOfClusters = numOfClusters, fName = saveFileName, savePath = plotBase)
}