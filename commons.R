
source("./R_script/util.R")
geneExpressionData = function(sourceFileName, sampleTypes, numOfClusters,states = 3 : 5,
	normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear"),
	discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile"),
	plotBase = "./plots",
	saveBase = "./processedData/discretizations") {

	set.seed(1024);
	plotPath = file.path(plotBase, sourceFileName);
	savePath = file.path(saveBase, sourceFileName);
	savePath = paste(c(savePath, "csv"), collapse = ".")

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

	for(i in 1 : length(discretization)) {
		titleName = discretization[i];
		bLine = baseLine(sourceData, numOfClusters = numOfClusters, normalization = match.fun(normalization[i]));
		similarity = pairWiseComparision(sourceData, numOfClusters = numOfClusters, benchMark = benchmark,
			states = states,
			normalization = normalization[i], discretization = match.fun(discretization[i]), dendro = dendro,
			savePath = plotPath,  titleName = paste(c(titleName,"dendrograph"), collapse = "_") )
		colnames(similarity) = paste(discretization[i], colnames(similarity), sep = " ")
		plotEvaluations(similarity, states = states, baseLine = bLine,
		 	titleName = titleName, savePath = plotPath, lx = 4.2, ly = 1)
		if(i == 1) {
			savedData = data.frame(similarity);
		}
		else {
			savedData = data.frame(savedData, similarity);
		}
	}
	write.csv(savedData, file = savePath)

}