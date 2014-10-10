
source("./R_script/util.R")


geneExpressionData = function(sourceFileName, sampleTypes, numOfClusters,states = 3 : 5,
	normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear"),
	discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile"),
	plotBase = "./plots",
	saveBase = "./processedData/discretizations") {

	set.seed(1024);
	plotPath = file.path(plotBase, sourceFileName);
	savePath = file.path(saveBase, sourceFileName);
	#savePath = paste(c(savePath, "csv"), collapse = ".")

	sData = read.csv(sourceFileName, sep = "\t", header=TRUE);
	hInfo = colnames(sData);
	hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));
	sourceData = sData[,hInfo > 0 & hInfo <= length(sampleTypes)];
	sourceData = as.matrix(sourceData);
	if(length(sampleTypes) > 1) {
		benchmark = hInfo[hInfo > 0 & hInfo  <= length(sampleTypes)];
		dendro = 3
		colnames(sourceData) = benchmark
	}
	else {
		benchmark = c()
		dendro = F
	}

	plotBarGraph(sourceData, numOfClusters = numOfClusters, fName = sourceFileName, savePath = plotPath)

}