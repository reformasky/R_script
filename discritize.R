#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)
#set seed for kmeans.
set.seed(1024);
DEBUG = TRUE;

#read data into R
#Should be read from a configuration file eventually. 
sourceData = read.csv(file = "D:/thesis/data/TCGA_BreastCancer_Centered.xls", header = FALSE, sep = '\t', nrows = 18055,colClasses = "numeric");
sourceData = data.matrix(sourceData)

#normalize to [0,1] range, for each row(the same gene across different samples)
sourceData = apply(sourceData, MARGIN = 1, function(vec) {
	maxVec = max(vec);
	minVec = min(vec);
	(vec -minVec) / (maxVec - minVec) ;
});
#don't forget to transpose
sourceData = t(sourceData)

#stdevs need to be use multiple times, so calculate and save;
stdevs = apply(sourceData, MARGIN = 1, FUN = sd);

#only keep genes with large variances for futher considerations.
#returns an array of indexs corresponding to the selected genes.
filtering = function(src = sourceData, ratioKept = 0.2) {
	if(!exists("stdevs"))
		stdevs = apply(src, MARGIN = 1, FUN = sd);
	threathold = quantile(stdevs, 1 - ratioKept);
	which(stdevs > threathold);
}

selectedList = filtering();
filteredData = sourceData[selectedList,]

#given the selected gene-sample matrix as filteredData, fit from kmeans method,
#group is the group number, find the gene index(with respect to the src) 
#closest(Euclidian distance) to the center of the cluster.
findCentroid = function(src = filteredData, fit, group) {

	selected = which(fit$cluster == group);

	if(DEBUG) {
		print(paste(c("lengthOfSelected= ",length(selected)), collapse = "" ));
		print(paste(c("groupNumber = ",group), collapse = "" ));
	}

	if(length(selected) > 1) {
		center = apply(src[selected,], MARGIN = 2, mean);
		distances = apply(src[selected,], MARGIN = 1, function(vec) {
			sum((vec - center) **2);
		}
		);
		selected[which.min(distances)];
	}
	else{
		selected[1];
	}

}

#selectGeneByKmeansCentriods, returns a list of indexs(with respect to srouceData, instead of filterdData),
# which coresponding to the 
#centroids of each cluster from k-means method.
selectGeneByKmeansCentriods = 
	function(src = filteredData, mapping = selectedList, k = 100) {
	fit = kmeans(src, k);
	result = lapply(1 : k, function(g) {
		findCentroid(src, fit, g);
		}
	);
	result = unlist(result);
	result = mapping[result];
}


#discritize
discritize = function(src, numOfLevels) {
	func = function(vec) {
		floor(vec * numOfLevels);
	}
	result = apply(src, MARGIN  = 1 ,func);
	t(result);
}

myPlot = function(folder = "D:/thesis/plots/kmeans", numOfGenes, numOfGroups,  x, y) {
	fileName = paste(c("numOfGroups", numOfGroups, "numOfGenes", numOfGenes), collapse = "_");
	savePath = paste(c(folder, fileName), collapse = "/");
	savePath = paste(c(savePath, "tiff"), collapse = ".");
	tiff(savePath);
	plot(x, y, main = fileName, xlab = "numOfLevels", ylab = "adjusted Rand index");
	dev.off();
}

#numOfClusters: the groups of samples for cutreee;
evaluation = function(numOfClusters) {
	numOfGenes = (1: 10) * 10;
	numOfLevels = (2 : 20) * 2;
	for(ng in numOfGenes) {
		selectedGenes = selectGeneByKmeansCentriods(k = ng);
		d = sourceData[selectedGenes,]
		distances = dist( t(d) );
		fit = hclust(distances, method = "ward");
		benchmark = cutree(fit, numOfClusters);
		tempResult = 1 : length(numOfLevels);
		for(i in  1: length(numOfLevels)) {
			distances = dist( t (discritize(d, numOfLevels[i])));
			fit = hclust(distances, method = "ward");
			cluster = cutree(fit, numOfClusters);
			tempResult[i] = RRand(cluster, benchmark) $adjRand;
			if(DEBUG) {
				print(paste(c("ng= ", ng, " numOfLevels= ", numOfLevels[i]), collapse = ""))
			}
		}
		myPlot(numOfGroups = numOfClusters, numOfGenes = ng, x = numOfLevels, y = tempResult);
		
	}
}

evaluation(10);

# #select genes by first filtering out genes with smaller variation; 
# #then select genes with smallest similarity with the one with largest variation
# #mode = {relative, absolute}. Relative mode filters out genes with lowest
# #threathold stdevs, and absolute mode filters out genes with stdevs lower than threathold 
# #data has been normalized into [0,1]
# #similarity is calculated basing on cosine similarity
# selectGeneBySimilarity = function(src = sourceData, threathold = 0.4, mode = "relative", numOfGenesKept = 500) {

# 	filteredGenes = stdevs > threathold;
# 	similarity = vector(mode = "numeric",length = length(stdevs));
# 	maxIndex = which.max(stdevs);
# 	maxLength = sqrt(src[maxIndex,] %*% src[maxIndex,]);
# 	for(i in 1 : length(stdevs)) {
# 		if(filteredGenes[i] == TRUE) {
# 			length = sqrt(src[i,] %*% src[i,]);
# 			similarity[i] = (src[i,] %*% src[maxIndex,]) / (maxLength * length)
# 		}
# 		else similarity[i] = 1;
# 	}
# 	similarity[maxIndex] = 0;
# 	order(similarity)[1 : numOfGenesKept];
# }




# evaluation = function(src, k = 10) {
# 	distances = dist(t(src));
# 	fit = hclust(distances, method = "ward");
# 	originalGroups = cutree(fit, k);
# 	steps = 1: 10;
# 	results = vector(mode = "numeric",length = length(steps));
# 	for(i in 1 : length(steps)) {
# 		distances = dist(t(discritize(src, steps[i] * 10)));
# 		fit = hclust(distances, method = "ward");
# 		groups = cutree(fit, k);
# 		results[i] = RRand(groups, originalGroups)$adjRand;
# 		print(results[i]);
# 	}
# 	results;
# }

# selectedVariance = selectGeneByVariation();
# selectedSimilarity = selectGeneBySimilarity();

# result = evaluation(sourceData[selectedSimilarity,]);

#plot the histogram of selected genes;

# plotSelected = function(mtx, prefix, folder = "D:/thesis/plots/") {
# 	savePath  = c(folder, prefix, "/","",".tiff")
# 	for( i in 1 : dim(mtx)[1]) {
# 		savePath[4] = i;
# 		name = paste(savePath, collapse = "");
# 		tiff(name)
# 		hist(mtx[i,], main = paste(c(prefix,"_", i), collapse = ""))
# 		dev.off()
# 	}
# }

# plotSelected(sourceData[selectedVariance,], prefix = "variation");
# plotSelected(sourceData[selectedSimilarity,], prefix = "similarity");