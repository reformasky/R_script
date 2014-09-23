
source("./R_script/util.R")
geneExpressionData = function(sourceFileName, sampleTypes, numOfClusters,states = 3 : 5,
	normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear"),
	discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile"),
	plotBase = "./plots") {

	set.seed(1024);
	basePath = file.path(plotBase, sourceFileName);

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
			savePath = basePath,  titleName = paste(c(titleName,"dendrograph"), collapse = "_") )
		plotEvaluations(similarity, states = states, baseLine = bLine,
		 	titleName = titleName, savePath = basePath, lx = 4.2, ly = 1)
	}

}