#extract the sample information from header file: for non-samples(i.e.,
# GeneId, and other statistic results) return -1; for pValues, return 0;
processHeader = function(str) {
	sampleTypes = c("Basal","Her2","LumA","LumB","Normal", "NA")
	stringVec = unlist(strsplit(str, split = ".", fixed = TRUE));
	if (length(stringVec) >= 2){
		which(sampleTypes == stringVec[1])
	}
	else if (str == "Pvalues"){
		0;
	}
	else {
		-1
	}
}

#linearly scale a vector linearly
normalizationLinear = function(vec) {
	minVec = min(vec);
	maxVec = max(vec);
	(vec - minVec)/ (maxVec - minVec);
}

discretizationFloor = function(vec, numOfStates = 3) {
	vec = normalizationLinear(vec);
	floor(vec * numOfStates);
} 


#normalize a vec basing on zScore; remember, when used in apply, if MARGIN == 1, will resulted a transposed matrix
normalizationZScore = function(vec) {
	sdVec = sd(vec);
	meanVec = mean(vec);
	(vec - meanVec) / sdVec;
}

#discretize a vec into numOfStates, wach state with equal probility
discretizationZScore = function(vec, numOfStates = 3) {
	floor( pnorm(normalizationZScore(vec)) * numOfStates );
}

#calculate the hamming distance between ROWs of matrix.
hammingMatrix = function(mtx) {
	hammingDistance = function(vec1, vec2) {
		sum(as.integer(vec1) != as.integer(vec2))
	}

	result = matrix(0, nrow = dim(mtx)[1], ncol = dim(mtx)[1])

	for( i in 1 : (dim(mtx)[1] - 1) ){
		for(j in (i + 1) : dim(mtx)[1]) {
			result[i,j] = result[j,i] = hammingDistance(mtx[i,], mtx[j,])
		}
	}
	as.dist(result)
}


baseLine = function(sData, numOfClusters, normalization = normalizationZScore) {
	sourceData = apply(sData, MARGIN = 1, normalization);
	distances = dist(sourceData)
	fit = hclust(distances)
	bMark = cutree(fit, k = numOfClusters)
	
	shuffleData = apply(t(sourceData), MARGIN = 1, sample)
	distances = dist(shuffleData)
	fit = hclust(distances)
	randomMark = cutree(fit, k = numOfClusters);
	baseLine = RRand(bMark, randomMark)$adjRand
}
# evaluate the effect of discretization on hclust, compares the clustering before and after discretization;
# if the true label is supplied, also compares the clustering of true label and after discretization;
# use two distance functions, euclidian and hamming distance.
pairWiseComparision = function(sData, numOfClusters, benchMark = c(),states = 3 : 10, discretization = discretizationZScore) {
	sourceData = apply(sData, MARGIN = 1, normalizationZScore);
	distances = dist(sourceData)
	fit = hclust(distances)
	bMark = cutree(fit, k = numOfClusters)
	
	cluster = function(sData, numOfStates) {
		euclideanGroups = sample(length(bMark), x = 1 : numOfClusters, replace = TRUE);
		hammingGroups = sample(length(bMark), x = 1 : numOfClusters, replace = TRUE);
		sourceData = apply(sData, MARGIN = 1, discretization, numOfStates);
		euclideanDistance = dist(sourceData);
		hammingDistance = hammingMatrix(sourceData);
		euclideanFit = hclust(euclideanDistance);
		hammingFit = hclust(hammingDistance);
		euclideanGroups = cutree(euclideanFit, k = numOfClusters);
		hammingGroups = cutree(hammingFit, k = numOfClusters);
		data.frame(euclideanGroups, hammingGroups);
	}

	euclideanPairWise = replicate(length(states), 0);
	hammingPairWise  = replicate(length(states), 0);
	euclideanBench = replicate(length(states), 0);
	hammingBench = replicate(length(states), 0)
	cat(paste(replicate(length(states), "*"), collapse = ""))
	print("")
	for(i in 1 : length(states)) {
		euclideanGroups = cluster(sData, states[i])$euclideanGroups;
		euclideanPairWise[i] = RRand(bMark, euclideanGroups)$adjRand;
		hammingGroups = cluster(sData, states[i])$hammingGroups;
		hammingPairWise[i] = RRand(bMark, hammingGroups)$adjRand;
		if(length(benchMark) != 0) {
			euclideanBench[i] = RRand(benchMark, euclideanGroups)$adjRand;
			hammingBench[i] = RRand(benchMark, hammingGroups)$adjRand;
		}
		cat(sprintf("-"))
	}
	print("");
	if(length(benchMark) != 0) {
		return(data.frame(euclideanPairWise, hammingPairWise, euclideanBench,hammingBench));
	}
	else {
		return(data.frame(euclideanPairWise, hammingPairWise));
	}
}

plotEvaluations = function(resultSet, states, baseLine, titleName, savePath, 
	xLab = "Number of States", yLab = "adjusted Rand index", lx = 3, ly = 0.2) {

	fileName = paste(c(titleName, "tiff"),collapse = ".")
	tiff(file.path(savePath, fileName));
	resultNumber = length(colnames(resultSet));
	for(i in 1 : resultNumber) {
		if( i == 1) {
			plot(states, unlist(resultSet[i]), xlab = xLab, ylab = yLab, 
				main = titleName,  ylim = c(baseLine - 0.05, 1), pch = i, cex = 1.5)
		}
		else {
			par(new = T)
			plot(states,unlist(resultSet[i]), xlab = xLab, ylab = yLab, ylim = c(baseLine - 0.05, 1), 
		xaxt = "n", yaxt="n", pch = i, cex = 1.5);
		}
	}
	abline(baseLine, 0, col = "blue", lty = 3);
	text(length(states)/2 + 2.5, baseLine + 0.05, "baseLine",col = "blue");
	legend(lx, ly, colnames(resultSet), pch=1 : resultNumber,
	 cex = replicate(resultNumber, 1), pt.cex = replicate(resultNumber, 1.5));
	dev.off(); 
	
}



# discretizationEqual = function(vec, numOfStates = 5) {
# 	step = ceiling(length(vec) / numOfStates);
#  	ranks = rank(vec);
#  	(ranks) %/% step
#  }


# evaluationFunctions = function(sData, numOfClusters, scalingMethod, distanceMethod, clusterMethod, benchMark) {
#  	groups = sample(length(benchMark), x = 1: numOfClusters, replace = TRUE);
#  	tryCatch({
#  		sourceData = apply(sData, MARGIN = 1, match.fun(scalingMethod));
# 		distances = dist(sourceData, method = distanceMethod);
# 		fit = hclust(distances, method = clusterMethod);
# 		groups = cutree(fit, k= numOfClusters);
#  	},
#  	error = function(cond){print("err")})
 	
# 	RRand(groups, benchMark);
# }

# plotFunctions = function(sourceData, savePath) {
# 	#parameters for hclust;
# 	scalingMethods = c("zeroOneScaling", "zscoreScaling");
# 	distanceMethods = c( "euclidean", "maximum", "manhattan",  "binary", "minkowski","canberra");
# 	clusterMethods = c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid");
# 	clusterLabels = c("wa","si","co","av","mc",",me","ce")
# 	labels = replicate( length(clusterMethods), "")
# 	results = replicate( length(clusterMethods), 1)

# 	for(s in 1 : length(scalingMethods)){
# 		for(d in 1 : length(distanceMethods)){
# 			cat(paste(replicate(length(clusterMethods), "*"), collapse = ""))
# 			print("")
# 			title = paste(c(scalingMethods[s], distanceMethods[d]), collapse = "_")
# 			counter = 1;
# 			for(cl in 1 : length(clusterMethods)) {
# 				label = clusterLabels[cl];
# 				result = evaluationFunctions(sourceData,  5,
# 				 	scalingMethods[s], distanceMethods[d],
# 				 	clusterMethods[cl], benchmark)$Rand;
# 				labels[counter] = label;
# 				results[counter] = result;
# 				counter = counter + 1;
# 				cat(sprintf("-"))
# 			}
# 			print("")
# 			fileName = paste(c(savePath, title), collapse = "/")
# 			fileName = paste(c(fileName, "tiff"), collapse = ".")
# 			tiff(fileName);
# 			plot(factor(labels), results, main= title, xlab= "agglomeration methods", ylab = "Rand index", ylim= c(0.2,0.8))
# 			dev.off();
# 		}		
# 	}
# }

# #test different number of states for zeroOneScaling
# plotStates = function(sourceData, savePath, numOfClusters, benchmark, range = 3 : 20) {
# 	evaluation = function(sData, numOfStates,  benchMark, numOfClusters,  distanceMethod = "euclidean", clusterMethod = "ward") {
#  		groups = sample(length(benchMark), x = 1: numOfClusters, replace = TRUE);
#  		tryCatch({
#  			sourceData = apply(sData, MARGIN = 1, discretizationFloor, numOfStates);
# 			distances = dist(sourceData, method = distanceMethod);
# 			fit = hclust(distances, method = clusterMethod);
# 			groups = cutree(fit, k= numOfClusters);
#  		},
#  		error = function(cond){print("err")})
# 		RRand(groups, benchMark);
# 	}

# 	titleName = paste(c("euclidean", "ward"), collapse = "_")
# 	fileName = paste(c(titleName, "tiff"),collapse = ".")
# 	results = replicate(length(range), 0);
# 	cat(paste(replicate(length(range), "*"), collapse = ""))
# 	print("")
# 	for(i in 1 : length(range)) {
# 		results[i] = evaluation(sourceData, range[i], benchmark, numOfClusters)$Rand;
# 		cat(sprintf("-"))
# 	}
# 	print("");
# 	tiff(paste(c(savePath,fileName), collapse = "/"));
# 	plot(range, results, xlab= "Number of States", ylab = "Rand index", main = titleName);
# 	dev.off();
# }