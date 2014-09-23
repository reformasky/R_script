library(phyclust, quiet = TRUE)

library(ggdendro)
library(ggplot2)
library(gridExtra)
#extract the sample information from header file (currently only good for genome portal data): for non-samples(i.e.,
# GeneId, and other statistic results) return -1; for pValues, return 0;
processHeader = function(str, sampleTypes) {
	
	stringVec = unlist(strsplit(str, split = ".", fixed = TRUE));
	if (length(stringVec) >= 2){
		which(sampleTypes == stringVec[1])
	}
	#only useful for selecting the features with small pValue
	else if (str == "Pvalues"){
		0;
	}
	# mark as discard data
	else {
		-1
	}
}

#linearly scale a vector linearly into [0,1]
normalizationLinear = function(vec) {
	minVec = min(vec);
	maxVec = max(vec);
	(vec - minVec)/ (maxVec - minVec);
}

discretizationFloor = function(vec, numOfStates = 3) {
	vec = normalizationLinear(vec);
	vec = floor(vec * numOfStates);
	vec[vec == numOfStates] = numOfStates - 1;
	vec;
} 


discretizationQuantile = function(vec, numOfStates = 3) {
	vec = rank(vec);
	vec = floor(vec * numOfStates / length(vec));
	vec[vec == numOfStates] = numOfStates - 1;
	vec;
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

#calculate the hamming distance between ROWs of matrix, returns a dist object
hammingMatrix = function(mtx) {
	hammingDistance = function(vec1, vec2) {
		sum(as.integer(vec1) != as.integer(vec2))
	}

	result = matrix(0, nrow = dim(mtx)[1], ncol = dim(mtx)[1])
	colnames(result) = rownames(result) = rownames(mtx);
	for( i in 1 : (dim(mtx)[1] - 1) ){
		for(j in (i + 1) : dim(mtx)[1]) {
			result[i,j] = result[j,i] = hammingDistance(mtx[i,], mtx[j,])
		}
	}
	as.dist(result)
}


# generate a baseLine by calculating the cluster assignment with orginal data and shuffled data(per gene)
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
# return a dataFrame 
pairWiseComparision = function(sData, numOfClusters, benchMark = c(),states = 3 : 10, normalization, 
	discretization, dendro = FALSE, savePath, titleName) {
	sourceData = apply(sData, MARGIN = 1, normalization);
	distances = dist(sourceData)
	fit = hclust(distances)
	if(dendro != FALSE) {
		plotDendrogram(fit, savePath = savePath, titleName = titleName);
	}

	bMark = cutree(fit, k = numOfClusters)
	
	euclidianCluster = function(sData, numOfStates) {
		sourceData = apply(sData, MARGIN = 1, discretization, numOfStates);
		euclideanDistance = dist(sourceData);
		euclideanFit = hclust(euclideanDistance);
	}

	hammingCluster = function(sData, numOfStates) {
		sourceData = apply(sData, MARGIN = 1, discretization, numOfStates);
		hammingDistance = hammingMatrix(sourceData);
		hammingFit = hclust(hammingDistance);
	}

	euclideanPairWise = replicate(length(states), 0);
	hammingPairWise  = replicate(length(states), 0);
	euclideanBench = replicate(length(states), 0);
	hammingBench = replicate(length(states), 0)
	cat(paste(replicate(length(states), "*"), collapse = ""))
	print("")
	for(i in 1 : length(states)) {
		euclideanFit = euclidianCluster(sData, states[i]);
		euclideanGroups = cutree(euclideanFit, k = numOfClusters);
		euclideanPairWise[i] = RRand(bMark, euclideanGroups)$adjRand;

		hammingFit = hammingCluster(sData, states[i])
		hammingGroups = cutree(hammingFit, k = numOfClusters);
		hammingPairWise[i] = RRand(bMark, hammingGroups)$adjRand;
		if(length(benchMark) > 1) {
			euclideanBench[i] = RRand(benchMark, euclideanGroups)$adjRand;
			hammingBench[i] = RRand(benchMark, hammingGroups)$adjRand;
		}

		if(dendro == states[i]) {
			plotDendrogram(euclideanFit, savePath = savePath, 
				titleName = paste(c(titleName, "euclidean", "numOfStates", dendro), collapse = "_"));
			plotDendrogram(hammingFit, savePath = savePath,
				titleName = paste(c(titleName, "hamming", "numOfStates", dendro), collapse = "_"))
		}
		cat(sprintf("-"))
	}
	print("");
	if(length(benchMark) > 1) {
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

#fit: from hclust; benchMark, true labeling
plotDendrogram = function(fit, savePath, titleName) {
	p1 = ggdendrogram(fit, rotate = FALSE, leaf_labels = FALSE, labels = FALSE);
	df = data.frame(seq = 1 : length(fit$labels), cluster = fit$labels[fit$order]);
	p2<-ggplot(df,aes(seq,y=1,fill=factor(cluster)))+geom_tile()+
	  scale_y_continuous(expand=c(0,0))+
	  theme(axis.title=element_blank(),
	        axis.ticks=element_blank(),
	        axis.text=element_blank(),
	        plot.background = element_blank(),
   			panel.grid.major = element_blank(),
   			panel.grid.minor = element_blank(),
   			panel.border = element_blank(),
	        legend.position="none");
	gp1<-ggplotGrob(p1)
	gp2<-ggplotGrob(p2)
	maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
	gp1$widths[2:5] <- as.list(maxWidth)
	gp2$widths[2:5] <- as.list(maxWidth)
	tiff(file.path(savePath, paste(c(titleName,".tiff"), collapse = "")));
	grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5), main = titleName); 
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