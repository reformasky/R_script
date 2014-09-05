#extract the sample information from header file: for non-samples(i.e.,
# GeneId, and other statistic results) return -1;
processHeader = function(str) {
	sampleTypes = c("Basal","Her2","LumA","LumB","Normal", "NA")
	stringVec = unlist(strsplit(str, split = ".", fixed = TRUE));
	if (length(stringVec) >= 2){
		which(sampleTypes == stringVec[1])
	}
	else{
		-1;
	}
}

zeroOneScaling = function(vec) {
	minVec = min(vec);
	maxVec = max(vec);
	(vec - minVec)/ (maxVec - minVec);
}

zscoreScaling = function(vec) {
	sdVec = sd(vec);
	meanVec = mean(vec);
	round((vec - meanVec) / sdVec);
}


evaluation = function(sData, numOfClusters, scalingMethod, distanceMethod, clusterMethod, benchMark) {
 	groups = benchMark;
 	sourceData = apply(sData, MARGIN = 1, match.fun(scalingMethod));
	distances = dist(sourceData, method = distanceMethod);
	fit = hclust(distances, method = clusterMethod);
	groups = cutree(fit, k= numOfClusters);

	RRand(groups, benchMark);
}