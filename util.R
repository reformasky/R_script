library(phyclust, quiet = TRUE)

library(ggdendro)
library(ggplot2)
library(gridExtra)
library(gplots)
#extract the sample information from header file (currently only good for genome portal data): for non-samples(i.e.,
# GeneId, and other statistic results) return -1; for pValues, return 0;
processHeader = function(str, sampleTypes) {
	
	stringVec = unlist(strsplit(str, split = ".", fixed = TRUE));
	if (length(stringVec) >= 2){
		if(stringVec[1] %in% sampleTypes)
				which(sampleTypes == stringVec[1])
		else
			-1
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
	vec = as.numeric(vec)
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
	vec = as.numeric(vec)
	sdVec = sd(vec);
	meanVec = mean(vec);
	(vec - meanVec) / sdVec;

}

#discretize a vec into numOfStates, wach state with equal probility
discretizationZScore = function(vec, numOfStates = 3) {
	result = floor( pnorm(normalizationZScore(vec)) * numOfStates );
	result[result >= numOfStates] = numOfStates - 1;
	result;
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


cosineMatrix = function(mtx) {
	cosineSimilarity = function(vec1 , vec2) {
		(vec1 %*% vec2) / sqrt(vec1 %*% vec1) / sqrt(vec2 %*% vec2)
	}
	result = matrix(1, nrow = dim(mtx)[1], ncol = dim(mtx)[1])
	colnames(result) = rownames(result) = rownames(mtx);
	for( i in 1 : (dim(mtx)[1] - 1) ){
		for(j in (i + 1) : dim(mtx)[1]) {
			result[i,j] = result[j,i] =1/ cosineSimilarity(mtx[i,], mtx[j,])
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
	discretization, dendro = FALSE, savePath = "./plots", titleName = "") {

	sourceData = apply(sData, MARGIN = 1, normalization);
	rownames(sourceData) = colnames(sData)
	distances = dist(sourceData)
	fit = hclust(distances)
	
	if(dendro != FALSE) {
		plotDendrogram(fit, savePath = savePath, titleName = titleName);
	}

	bMark = cutree(fit, k = numOfClusters)
	
	euclidianCluster = function(sData, numOfStates) {
		sourceData = apply(sData, MARGIN = 1, discretization, numOfStates);
		rownames(sourceData) = colnames(sData)

		euclideanDistance = dist(sourceData);
		euclideanFit = hclust(euclideanDistance);
	}

	hammingCluster = function(sData, numOfStates) {
		sourceData = apply(sData, MARGIN = 1, discretization, numOfStates);
		rownames(sourceData) = colnames(sData)

		hammingDistance = hammingMatrix(sourceData);
		hammingFit = hclust(hammingDistance);
	}

	eVsNoDiscretize = replicate(length(states), 0);
	hvsNoDiscretize  = replicate(length(states), 0);
	evsTrueLabel = replicate(length(states), 0);
	hvsTrueLabel = replicate(length(states), 0)
	cat(paste(replicate(length(states), "*"), collapse = ""))
	print("")
	for(i in 1 : length(states)) {
		euclideanFit = euclidianCluster(sData, states[i]);
		euclideanGroups = cutree(euclideanFit, k = numOfClusters);
		eVsNoDiscretize[i] = RRand(bMark, euclideanGroups)$adjRand;

		hammingFit = hammingCluster(sData, states[i])
		hammingGroups = cutree(hammingFit, k = numOfClusters);
		hvsNoDiscretize[i] = RRand(bMark, hammingGroups)$adjRand;
		if(length(benchMark) > 1) {
			evsTrueLabel[i] = RRand(benchMark, euclideanGroups)$adjRand;
			hvsTrueLabel[i] = RRand(benchMark, hammingGroups)$adjRand;
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
		discretized  = data.frame(eVsNoDiscretize, hvsNoDiscretize, evsTrueLabel,hvsTrueLabel);

	}
	else {
		discretized = data.frame(eVsNoDiscretize, hvsNoDiscretize);
	}
	rownames(discretized) = states;
	discretized = t(discretized)
	noDiscretizeVsTrueLabel = 1;
	if(length(benchMark) > 0)
		noDiscretizeVsTrueLabel = RRand(benchMark, bMark) $adjRand
	result = list(discretized = discretized, noDiscretizeVsTrueLabel = noDiscretizeVsTrueLabel )
	return(result);
}

plotBarGraph = function(sourceData, numOfClusters, fName, savePath,
	normalization = c("normalizationLinear", "normalizationLinear", "normalizationZScore"), 
 	discretization = c( "discretizationQuantile",  "discretizationFloor", "discretizationZScore"), 
 	states = 3 : 5) {
	basePath = "./plots"
	benchMark = as.numeric(colnames(sourceData))

	tiff(file.path(savePath, paste(c(fName, ".tiff"), collapse = "")), units="in", width=8, height=6, res=300)
	par(mfcol  = c(2, length(discretization)), mar = c(2,5,1,1))


	for(i in 1 : length(discretization)) {
		titleName = paste(c(discretization[i],fName), collapse = "_")
		# bLine = baseLine(sourceData, numOfClusters, normalization = match.fun(normalization[i]))
		results = pairWiseComparision(sourceData, numOfClusters = numOfClusters, states = states, 
			normalization = normalization[i], discretization = match.fun(discretization[i]), benchMark = benchMark, 
			dendro = F, savePath = basePath,titleName = paste(c(titleName,"dendrograph"), collapse = "_"));
		referenceLine = unlist(results$noDiscretizeVsTrueLabel)[1]
		discretized = unlist(results$discretized)
		plotData = cbind(discretized[1:2], discretized[5:6], discretized[9:10])
		if(i == 1) {
			bp = barplot(plotData,beside = T, col = c("gray", "black"), xaxt = "n",  ylim = c(-0.05,1), ylab = "adjusted Rand Index")
		}
		else {
			bp = barplot(plotData,beside = T, col = c("gray", "black"), xaxt = "n", yaxt = "n",  ylim = c(-0.05,1))
		}		

		plotData = cbind(discretized[3:4], discretized[7:8], discretized[11:12])
		if(i == 1) {
			bp = barplot(plotData,beside = T, col = c("red", "blue"), xaxt = "n",  ylim = c(-0.05,1), ylab = "adjusted Rand Index")
		}
		else {
			bp = barplot(plotData,beside = T, col = c("red", "blue"), xaxt = "n", yaxt = "n",  ylim = c(-0.05,1))
		}
		axis(1, at= c( 2, 5, 8), labels = c( 3, 4,5), tick = F)

		abline(referenceLine, 0, col = "green", lty = 3, lwd = 2)
		if(i == 1) {
			saved = c(results)
		}
		else{
			saved = c(saved, results)
		}
	}
	dev.off()
	saved
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
	tiff(file.path(savePath, paste(c(titleName,".tiff"), collapse = "")), units="in", width=11, height= 8, res=300);
	grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5), main = titleName); 
	dev.off(); 
}




evaluationFunctions = function(sData, numOfClusters, scalingMethod, distanceMethod, clusterMethod) {
	benchMark = as.numeric(colnames(sData))
 	groups = sample(length(benchMark), x = 1: numOfClusters, replace = TRUE);
 	tryCatch({
 		sourceData = apply(sData, MARGIN = 1, match.fun(scalingMethod));
		distances = dist(sourceData, method = distanceMethod);
		fit = hclust(distances, method = clusterMethod);
		groups = cutree(fit, k= numOfClusters);
 	},
 	error = function(cond){print(scalingMethod);
 		print(distanceMethod);
 		print(clusterMethod)})
 	
	RRand(groups, benchMark);
}

evaluateDistanceAndClustering = function(sData, numOfClusters = 4, titleName, scalingMethod,
	distances = c( "euclidean", "manhattan"),
	clusterings = c( "single", "complete", "average", "median", "centroid", "ward"),
	savePath = "./plots/distancesAndClusterings") {
	xLab = "Clustering Algorithm"
	yLab = "adjusted Rand Index"
	savePath = file.path(savePath, paste(c(titleName, "tiff"), collapse = "."));
	result = matrix(0, nrow = length(distances), ncol = length(clusterings));
	rownames(result) = distances;
	colnames(result) = clusterings;
	benchmark = as.numeric(rownames(sData));
	for(r in 1 : dim(result)[1]) {
		for(c in 1 : dim(result)[2]) {
			
			result[r,c] = evaluationFunctions(sData, numOfClusters = numOfClusters, scalingMethod = scalingMethod,
				distanceMethod = distances[r], clusterMethod = clusterings[c])$adjRand;
		}
	}
	# result = format(result, digits = 3)
	# tiff(file=savePath, units="in", width=5.5, height= 4, res=300);
	# plot(1 : dim(result)[2], result[1,], xlab = xLab, ylab = yLab, xaxt = "n",
	# 			main = titleName,  ylim = c(-0.05, 1), pch = 1, cex = 1.5)
	# for(i in 2 : dim(result)[1]) {
	# 	par(new = T)
	# 	plot(1 : dim(result)[2],result[i,], xlab = xLab, ylab = yLab, ylim = c( - 0.05, 1), 
	# 	xaxt = "n", yaxt="n", pch = i, cex = 1.5)
	# }
	# axis(1, at = 1 : dim(result)[2], labels = colnames(result))
	# dev.off();
	result;
}

# discretizationEqual = function(vec, numOfStates = 5) {
# 	step = ceiling(length(vec) / numOfStates);
#  	ranks = rank(vec);
#  	(ranks) %/% step
#  }




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



# plotEvaluations = function(results, states, baseLine, titleName, savePath, 
# 	xLab = "Number of States", yLab = "adjusted Rand index", lx = 3, ly = 0.2) {

# 	resultSet = data.frame(results$discretized);
# 	noDiscretizeVsTrueLabel = unlist(results$noDiscretizeVsTrueLabel)[1];

# 	fileName = paste(c(titleName, "tiff"),collapse = ".")
# 	tiff(file.path(savePath, fileName),units="in", width=6, height=6, res=300);
# 	resultNumber = length(colnames(resultSet));
# 	color = "black"
# 	for(i in 1 : resultNumber) {
# 		yData = unlist(resultSet[i])
# 		if( i == 1) {
# 			plot(states, yData, xlab = xLab, ylab = yLab, xaxt = "n",
# 				main = titleName,  ylim = c(baseLine - 0.05, 1), pch = i, cex = 1.5, color = color)
# 		}
# 		else {
# 			if ( i > (resultNumber / 2 ))
# 			color = "red"
# 			par(new = T)
# 			plot(states,yData, xlab = xLab, ylab = yLab, ylim = c(baseLine - 0.05, 1), 
# 		xaxt = "n", yaxt="n", pch = i, cex = 1.5, col = color);
# 		}
# 	}
# 	axis(1, at = states)
# 	abline(baseLine, 0, col = "blue", lty = 3)
# 	if(noDiscretizeVsTrueLabel < 1) {
# 		abline(noDiscretizeVsTrueLabel, 0, col = "red", lty = 5);

# 	}

# 	# legend(lx, ly, colnames(resultSet), pch=1 : resultNumber,
# 	#  cex = replicate(resultNumber, 1), pt.cex = replicate(resultNumber, 1.5));
# 	dev.off(); 
	
# }