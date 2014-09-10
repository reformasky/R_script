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

# only support 5 states for now; the numOfStates is a fake argument 
#just to make it compatible with "discretizationFloor"
discretizationZscore = function(vec, numOfStates= 5) {
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

evaluationFunctions = function(sData, numOfClusters, scalingMethod, distanceMethod, clusterMethod, benchMark) {
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

plotFunctions = function(sourceData, savePath) {
	#parameters for hclust;
	scalingMethods = c("zeroOneScaling", "zscoreScaling");
	distanceMethods = c( "euclidean", "maximum", "manhattan",  "binary", "minkowski","canberra");
	clusterMethods = c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid");
	clusterLabels = c("wa","si","co","av","mc",",me","ce")
	labels = replicate( length(clusterMethods), "")
	results = replicate( length(clusterMethods), 1)

	for(s in 1 : length(scalingMethods)){
		for(d in 1 : length(distanceMethods)){
			cat(paste(replicate(length(clusterMethods), "*"), collapse = ""))
			print("")
			title = paste(c(scalingMethods[s], distanceMethods[d]), collapse = "_")
			counter = 1;
			for(cl in 1 : length(clusterMethods)) {
				label = clusterLabels[cl];
				result = evaluationFunctions(sourceData,  5,
				 	scalingMethods[s], distanceMethods[d],
				 	clusterMethods[cl], benchmark)$Rand;
				labels[counter] = label;
				results[counter] = result;
				counter = counter + 1;
				cat(sprintf("-"))
			}
			print("")
			fileName = paste(c(savePath, title), collapse = "/")
			fileName = paste(c(fileName, "tiff"), collapse = ".")
			tiff(fileName);
			plot(factor(labels), results, main= title, xlab= "agglomeration methods", ylab = "Rand index", ylim= c(0.2,0.8))
			dev.off();
		}		
	}
}

#test different number of states for zeroOneScaling
plotStates = function(sourceData, savePath, numOfClusters, benchmark, range = 3 : 20) {
	evaluation = function(sData, numOfStates,  benchMark, numOfClusters,  distanceMethod = "euclidean", clusterMethod = "ward") {
 		groups = sample(length(benchMark), x = 1: numOfClusters, replace = TRUE);
 		tryCatch({
 			sourceData = apply(sData, MARGIN = 1, discretizationFloor, numOfStates);
			distances = dist(sourceData, method = distanceMethod);
			fit = hclust(distances, method = clusterMethod);
			groups = cutree(fit, k= numOfClusters);
 		},
 		error = function(cond){print("err")})
		RRand(groups, benchMark);
	}

	titleName = paste(c("euclidean", "ward"), collapse = "_")
	fileName = paste(c(titleName, "tiff"),collapse = ".")
	results = replicate(length(range), 0);
	cat(paste(replicate(length(range), "*"), collapse = ""))
	print("")
	for(i in 1 : length(range)) {
		results[i] = evaluation(sourceData, range[i], benchmark, numOfClusters)$Rand;
		cat(sprintf("-"))
	}
	print("");
	tiff(paste(c(savePath,fileName), collapse = "/"));
	plot(range, results, xlab= "Number of States", ylab = "Rand index", main = titleName);
	dev.off();
}