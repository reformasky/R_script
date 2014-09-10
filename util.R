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

zeroOneScaling = function(vec, fun) {
	minVec = min(vec);
	maxVec = max(vec);
	(vec - minVec)/ (maxVec - minVec);
}

zscoreScaling = function(vec) {
	sdVec = sd(vec);
	meanVec = mean(vec);
	(vec - meanVec) / sdVec;
}

discretizationFloor = function(vec, numOfStates = 5) {
	vec = zeroOneScaling(vec);
	floor(vec * numOfStates);
} 

# only support 5 states for now
discretizationZscore = function(vec) {
	upperBound = c(-2, -0.5, 0.5, 2);
	vec = zscoreScaling(vec);
	filter = function(var) {
		for(i in 1 : length(upperBound)) {
			if (var < upperBound[i])
				return(i)
		}
		return(5)
	}
	unlist(lapply(vec, filter))
}


discretizationEqual = function(vec, numOfStates = 5) {
	step = ceiling(length(vec) / numOfStates);
	ranks = rank(vec);
	(ranks) %/% step
}

evaluation = function(sData, numOfClusters, scalingMethod, distanceMethod, clusterMethod, benchMark) {
 	groups = sample(length(benchMark), x = 1: numOfClusters, replace = TRUE);
 	tryCatch({
 		sourceData = apply(sData, MARGIN = 1, match.fun(scalingMethod));
		distances = dist(sourceData, method = distanceMethod);
		fit = hclust(distances, method = clusterMethod);
		groups = cutree(fit, k= numOfClusters);
 	},
 	error = function(cond){print("err")})
 	
	RRand(groups, benchMark);
}
