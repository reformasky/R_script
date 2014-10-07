
setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)
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
	benchMark = as.numeric(colnames(sourceData))
	plotBase = file.path(plotPath, paste(c("sd", sd), collapse = "_"))
	saveFileName =  paste(c("sd", sd), collapse = "_");
	
	
	for(i in 1 : length(discretization)) {
		titleName = paste(c(discretization[i],"simulation", "sd", sd), collapse = "_")
		bLine = baseLine(sourceData, numOfClusters, normalization = match.fun(normalization[i]))
		result = pairWiseComparision(sourceData, numOfClusters = numOfClusters, states = states, 
			normalization = normalization[i], discretization = match.fun(discretization[i]), benchMark = benchMark, 
			dendro = 3, savePath = plotBase,titleName = paste(c(titleName,"dendrograph"), collapse = "_"));
		plotEvaluations(result, states, baseLine = bLine, 
			titleName = titleName, savePath = plotBase, lx = 4.2, ly = 0.92)
		similarity = data.frame(result$discretized)
		noDiscretizeVsTrueLabel = result$noDiscretizeVsTrueLabel
		colnames(similarity) = paste(discretization[i], colnames(similarity), sep = " ")
		if(i == 1) {
			savedData = data.frame(similarity);
			noDiscretization = data.frame(noDiscretizeVsTrueLabel)
		}
		else {
			savedData = data.frame(savedData, similarity);
			noDiscretization = data.frame(noDiscretization, noDiscretizeVsTrueLabel )
		}
	}
	write.csv(noDiscretization, file = file.path(saveBase, paste(c("noD_",saveFileName, ".csv"), collapse = "")))
	saveFileName = paste(c(saveFileName, ".csv"), collapse = "")
	write.csv(savedData, file =  file.path(saveBase, saveFileName))

}