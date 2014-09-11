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


normalizationZScore = function(vec) {
	sdVec = sd(vec);
	meanVec = mean(vec);
	(vec - meanVec) / sdVec;
}

#discretize a vec into numOfStates
discretizationZScore = function(vec, numOfStates = 3) {
	floor( pnorm(normalizationZScore(vec)) * numOfStates);
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

pairWiseComparision = function(sData, savePath, numOfClusters, benchMark = c(),states = 3 : 10) {
	sourceData = apply(sData, MARGIN = 1, normalizationZScore);
	distances = dist(sourceData)
	fit = hclust(distances)
	bMark = cutree(fit, k = numOfClusters)

	shuffleData = apply(t(sourceData), MARGIN = 1, sample)
	distances = dist(shuffleData)
	fit = hclust(distances)
	randomMark = cutree(fit, k = numOfClusters);

	baseLine = RRand(bMark, randomMark)$adjRand

	cluster = function(sData, numOfStates) {
		euclideanGroups = sample(length(bMark), x = 1 : numOfClusters, replace = TRUE);
		hammingGroups = sample(length(bMark), x = 1 : numOfClusters, replace = TRUE);
		sourceData = apply(sData, MARGIN = 1, discretizationZScore, numOfStates);
		euclideanDistance = dist(sourceData);
		hammingDistance = hammingMatrix(sourceData);
		euclideanFit = hclust(euclideanDistance);
		hammingFit = hclust(hammingDistance);
		euclideanGroups = cutree(euclideanFit, k = numOfClusters);
		hammingGroups = cutree(hammingFit, k = numOfClusters);
		data.frame(euclideanGroups, hammingGroups);
	}

	euclideanCluster = function(sData, numOfStates) {
 		groups = sample(length(bMark), x = 1: numOfClusters, replace = TRUE);
 		tryCatch({
 			sourceData = apply(sData, MARGIN = 1, discretizationZScore, numOfStates);
			distances = dist(sourceData);
			fit = hclust(distances, method = "ward");
			groups = cutree(fit, k= numOfClusters);
 		},
 		error = function(cond){print("err")})
 		
		groups;
	}

	hammingCluster = function(sData, numOfStates) {
		groups = sample(length(bMark), x = 1 : numOfClusters, replace = TRUE);
		tryCatch({
				sourceData = apply(sData, MARGIN = 1, discretizationZScore, numOfStates);
				distances = hammingMatrix(sourceData);
				fit = hclust(distances, method = "ward");
				groups = cutree(fit, k = numOfClusters);
			},
			error = function(cond){print("err")})
		groups;

	}


	euclideanPairWise = replicate(length(states), 0);
	hammingPairWise  = replicate(length(states), 0);
	euclideanPam50 = replicate(length(states), 0);
	hammingPam50 = replicate(length(states), 0)
	cat(paste(replicate(length(states), "*"), collapse = ""))
	print("")
	for(i in 1 : length(states)) {
		# euclideanGroups = euclideanCluster(sData, states[i])
		# euclideanPairWise[i] = RRand(bMark, euclideanGroups)$adjRand;
		euclideanGroups = cluster(sData, states[i])$euclideanGroups;
		euclideanPairWise[i] = RRand(bMark, euclideanGroups)$adjRand;
		hammingGroups = cluster(sData, states[i])$hammingGroups;
		hammingPairWise[i] = RRand(bMark, hammingGroups)$adjRand;
		if(length(benchMark) == length(bMark)) {
			euclideanPam50[i] = RRand(benchMark, euclideanGroups)$adjRand;
			hammingPam50[i] = RRand(benchMark, hammingGroups)$adjRand;
		}
		# tiff(paste(c(savePath,"/",states[i],".tiff"), collapse = ""));
		# plot(table(bMark, euclideanGroups), 
		# 	main = paste(c("numOfStates= ", states[i])), xlab = "before discretization", ylab = "after discretization")
		# dev.off();
		cat(sprintf("-"))
	}
	print("");
	# titleName = paste(c("euclidean", "ward"), collapse = "_")
	# fileName = paste(c(titleName, "tiff"),collapse = ".")
	# tiff(paste(c(savePath,fileName), collapse = "/"));
	# plot(states, euclideanPairWise, xlab= "Number of States", ylab = "Rand index", main = titleName, ylim = c(baseLine - 0.05, 0.7));
	# abline(baseLine, 0);
	# text(length(states)/2, baseLine + 0.05, "baseLine",col = "blue")
	# dev.off();
	titleName = "pairwise comparion_zscore discretization"
	fileName = paste(c(titleName, "tiff"),collapse = ".")
	tiff(paste(c(savePath,fileName), collapse = "/"));
	plot(states, euclideanPairWise, xlab = "Number of States", ylab = "Adjusted Rand Index", 
		main = titleName,  ylim = c(baseLine - 0.05, 0.75), pch = 15, col = "red")
	par(new = T)
	plot(states, hammingPairWise, xlab = "Number of States", ylab = "Adjusted Rand Index", ylim = c(baseLine - 0.05, 0.75), 
		xaxt = "n", yaxt="n", pch = 15, col = "green");
	
	if(length(benchMark) == length(bMark)) {
		par(new = T);
		plot(states, euclideanPam50, xlab = "Number of States", ylab = "Adjusted Rand Index", ylim = c(baseLine - 0.05, 0.75), 
		xaxt = "n", yaxt="n", pch = 16, col = "red");
		par(new = T);
		plot(states, hammingPam50, xlab = "Number of States", ylab = "Adjusted Rand Index", ylim = c(baseLine - 0.05, 0.75), 
		xaxt = "n", yaxt="n", pch = 16, col = "green");
	}
	abline(baseLine, 0, col = "blue");
	text(length(states)/2 + 2.5, baseLine + 0.05, "baseLine",col = "blue");
	if(length(benchMark) == length(bMark)) {
		legend(3, 0.2, c("euclidean_pair", "hamming_pair", "euclidean_pam50", "hamming_pam50"), pch=c(15,15, 16,16), col = c("red", "green", "red", "green"))
	}
	else {
		legend(3, 0.15, c("euclidean_pair", "hamming_pair"), pch=c(15,15), col = c("red", "green"))
	}
	dev.off(); 
	if(length(benchMark) != length(bMark))
		data.frame(euclideanPairWise, hammingPairWise)
	else 
		data.frame(euclideanPairWise, hammingPairWise, euclideanPam50, hammingPam50)
}










# zeroOneScaling = function(vec) {
# 	minVec = min(vec);
# 	maxVec = max(vec);
# 	(vec - minVec)/ (maxVec - minVec);
# }



# discretizationFloor = function(vec, numOfStates = 5) {
# 	vec = zeroOneScaling(vec);
# 	floor(vec * numOfStates);
# } 



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