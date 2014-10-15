
source("./R_script/util.R")

# returns the list of lists of genes with gene lists, each group as a seperate lists
gpLists = function(fileName) {
  lines = readLines(file(fileName, open = "r"))
  results = list()
  for( i in 1 : length(lines)) {
    tempList = unlist(strsplit(lines[i], ","))
    result = tempList[ tempList != "" ]
    results = c(results, list(result))
  }
  results
}


# returns a matrix which is a subset of sData, containing all geneIds appear in lists
  # the reason doing this is that we have to calculate 1D-Jury for each gene entry, 
  # and reducing the number of genes will make the program run faster

findCombinedList = function(sData, geneLists, maxGenesKept = 20, normalization = function(vec) {vec}) {
  geneIds = rownames(sData)

  selectedGeneId = c()
  for(geneList in geneLists) {
    geneList = unlist(geneList)
    indexes = intersect(geneList, geneIds)

    #selected = filterGenes(sData[indexes,], normalization, maxGenesKept)
    selectedGeneId = union(selectedGeneId, indexes)
  }

  indexes = match(selectedGeneId, geneIds)
  sData[indexes,]

}

#if a list contains too many genes, we only select the gene with highest varience
filterGenes = function(sData, normalization = function(vec) {vec}, maxGenesKept = 20 ) {
  #TODO: better define selected genes
  if (dim(sData)[1] <= maxGenesKept)
    return(rownames(sData))
  else{
    normalizedData = t( apply(sData, MARGIN = 1, normalization))
    sds = apply(normalizedData, MARGIN = 1, sd)
    names(sort(sds, decreasing = T)[1 : maxGenesKept])
  }
}


# return a list of lists of indexes, each list coresponds to one list in lists
findIndexesFromCombined = function(sData, lists) {
  results = list();
  vec = rownames(sData)
  for( i in 1 : length(lists) ) {
    lhs = unlist(lists[i])
    temp = match(lhs, vec);
    result = temp[!is.na(temp)]
    if(length(result) > 0)
    	results = c(results, list(result))
  }
  results
}

# assuming mtx is a matrix with proper discretization; returns nStates x numOfTotalGenes;
# mtx is transposed because of discretization! numOfSamples x numOfGenes
oneDJuryCache = function(mtx, numOfStates) {
	result = matrix(0, ncol = dim(mtx)[2], nrow = numOfStates);
	rownames(result) = 1 : numOfStates ;
	colnames(result) = colnames(mtx);
	for(r in 1: dim(result)[1]) {
		for(c in 1 : dim(result)[2]) {
			result[r, c] = sum(mtx[,c] == r );
		}
	}
	result;
}


# sData: matrix of discretized expression ; samples X gene;
# genesInGroups(very important): a list of lists of genes, each inner list is the geneIds 
	# for that particular groups of genes; need to use unlist to retrive an array of genes
# returns the numOfSample x numOfGroups matrix, each row is a 1DJury score vector.
oneDJuryScore = function(sData, genesInGroups, numOfStates = 3) {

	numOfGroups = length(genesInGroups);
	sourceData = sData + 1
	cachedStates = oneDJuryCache(sourceData, numOfStates);

	#translate [state, gene] score -> [sample, gene] score
	cachedSamples = matrix(0, nrow = dim(sourceData)[1], ncol = dim(sourceData)[2])
	for(r in 1 : dim(cachedSamples)[1]){
		for( c in 1 : dim(cachedSamples)[2]){
			cachedSamples[r,c] = cachedStates[sourceData[r,c], c];
		}
	}
	# numOfSamples x numOfGroups
	result = matrix(0, nrow = dim(sourceData)[1], ncol = numOfGroups);
	rownames(result) = rownames(sourceData);
	colnames(result) = 1 : numOfGroups
	for(group in 1 : numOfGroups) {
		selectedGenes = unlist(genesInGroups[group])
		result[,group] = t(rowSums(cachedSamples[, selectedGenes]))/ length(selectedGenes)
	}
	
	result;
}


evaluateOneDJury = function(sData, genesInGroups, normalization,  discretization, numOfClusters,
	numOfStates = c(3) , dendro = F,  plotBase = F) {
	titleName = paste(c(normalization, "noDiscretization","dengro"), collapse = "_")
	trueLabeling = as.numeric(colnames(sData));
	sourceData = apply(sData, MARGIN = 1, match.fun(normalization));
	rownames(sourceData) = as.numeric(colnames(sData))
	distance = dist(sourceData);
	fit = hclust(distance)
	nonDiscretizedGroups = cutree(fit, k = numOfClusters)
	plotDendrogram(fit, savePath = plotBase, 
	 	titleName = titleName);

	euclideanVsTrueLable = replicate(length(numOfStates), 0);
	euclieanVsNoDiscretize = replicate(length(numOfStates), 0);


	for(i in 1 : length(numOfStates)) {
		state = numOfStates[i]
		sourceData = apply(sData, MARGIN = 1, discretization , numOfStates = state)
		sourceData = oneDJuryScore(sData = sourceData, genesInGroups = genesInGroups, numOfStates = state);
		rownames(sourceData) = as.numeric(colnames(sData))
		distance = dist(sourceData);
		fit = hclust(distance)
		if(state == dendro){
			plotDendrogram(fit, savePath = plotBase, 
				titleName = paste(c(discretization, "numOfStates", state), collapse = "_"))
		}
		groups = cutree(fit, k = numOfClusters);
		euclideanVsTrueLable[i] = RRand(trueLabeling, groups)$adjRand
		euclieanVsNoDiscretize[i] = RRand(nonDiscretizedGroups, groups)$adjRand;

	}
	result = data.frame( vsNoDiscretization = euclieanVsNoDiscretize, vsTrueLabel = euclideanVsTrueLable)
	rownames(result) = numOfStates;
	result = format(result, digits = 3);
	result = list(discretized = result, noDiscretizeVsTrueLabel = RRand(trueLabeling, nonDiscretizedGroups)$adjRand)
}


#Assuming that all resultDump follows this pattern:
	# against_trueLabel = data.frame(original = results_sourceData, quantile = results_Quantile, interval = results_Floor, zScore = results_zScore,
	# 					oneD_feq = results_oneD_q, oneD_int = results_oneD_i, oneD_z = results_oneD_z)

	# against_noOneD = data.frame(feq = RRand(groups_Quantile, groups_oneD_q)$adjRand,
	# 							interval = RRand(groups_Floor, groups_oneD_i)$adjRand,
	# 							zScore = RRand(groups_zScore, groups_oneD_z)$adjRand)
	# result = list(vsTrueLabel = against_trueLabel, against_noOneD = against_noOneD)


plotOneDJury = function(fileName, loading = T) {
	load(file.path("processedData/oneDJury", fileName))
	against_trueLabel = unlist(result$vsTrueLabel)
	against_noOneD = unlist(result$against_noOneD)

	tiff(file.path("plots/oneDJury/geneLists",paste(c(fileName, ".tiff"), collapse = "" ) ), units="in", width=7, height=8, res=300)
	plotData = cbind(eq_feq = c(against_trueLabel[2], against_trueLabel[5], against_noOneD[1]),
					 eq_int = c(against_trueLabel[3], against_trueLabel[6], against_noOneD[2]),
					 z_score =c(against_trueLabel[4], against_trueLabel[7], against_noOneD[3]))
	barplot(plotData, beside = T, col = c("red", "green", "blue"), main = fileName, xaxt = "n",  ylim = c(-0.1,1), ylab = "adjusted Rand Index")
	axis(1, at = c(2.5,6.5,10.5), tick = F, labels = c("Equal_frequency", "Equal_interval", "Z Score"))
	abline(against_trueLabel[1], 0, lty = 3, lwd = 2, col = "black")
	dev.off()
}


plotOneDJuryBar = function(fileName, dataBase = "./processedData/oneDJury", plotBase = "./plots/oneDJury") {

	sData = read.csv(file.path(dataBase, paste(c(fileName, "csv"), collapse = "." ) ), header = T, sep = "," )
	savedFile = file.path(plotBase, paste(c(fileName, "tiff"), collapse = "." ))

	tiff(savedFile, units="in", width=8, height=6, res=300)
	par(mfcol  = c(2,3), mar = c(2,5,1,1))
	for(d in 1 : 3) {
		base = -1
		plotData = t(as.matrix(sData[, base + (3 *d)]))
		referenceLine = sData[1, d *3 + 1]
		if(d == 1) {
			bp = barplot(plotData,beside = T, col = c("black"), xaxt = "n",  ylim = c(-0.1,1), ylab = "adjusted Rand Index")
		}
		else {
			bp = barplot(plotData,beside = T, col = c("black"), xaxt = "n", yaxt = "n",  ylim = c(-0.1,1))
		}

		base = 0
		plotData = t(as.matrix(sData[, base + (3 *d)]))
		if(d == 1) {
			bp = barplot(plotData,beside = T, col = c("red"), xaxt = "n",  ylim = c(-0.1,1), ylab = "adjusted Rand Index")
		}
		else {
			bp = barplot(plotData,beside = T, col = c("red"), xaxt = "n", yaxt = "n",  ylim = c(-0.1,1))
		}
		axis(1, at= c( 1.5, 3.5, 5.5), labels = c( 3, 4,5), tick = F)

		abline(referenceLine, 0, col = "blue", lty = 3, lwd = 2)
	}
	dev.off()
}

# }

# evaluateOneDJury = function(sourceData, trueLabel, geneIndexs, numOfClusters, numOfState, fileName){

# 	distances = dist(t(sourceData))
# 	fit_sourceData = hclust(distances)
# 	groups_sourceData = cutree(fit_sourceData, k = numOfClusters)
# 	results_sourceData = RRand(groups_sourceData, trueLabel)$adjRand

# 	distances = dist(apply(sourceData, MARGIN = 1, discretizationZScore, numOfStates))
# 	fit_zScore = hclust(distances)
# 	groups_zScore = cutree(fit_zScore, k = numOfClusters)
# 	results_zScore = RRand(groups_zScore, trueLabel)$adjRand

# 	distances = dist(apply(sourceData, MARGIN = 1, discretizationFloor, numOfStates))
# 	fit_Floor = hclust(distances)
# 	groups_Floor = cutree(fit_Floor, k = numOfClusters)
# 	results_Floor = RRand(groups_Floor, trueLabel)$adjRand

# 	distances = dist(apply(sourceData, MARGIN = 1, discretizationQuantile, numOfStates))
# 	fit_Quantile = hclust(distances)
# 	groups_Quantile = cutree(fit_Quantile, k = numOfClusters)
# 	results_Quantile = RRand(groups_Quantile, trueLabel)$adjRand

# 	#oneDJuryScores = oneDJuryScore(sourceData, geneIndexs)

# 	#results = evaluateOneDJury(sourceData, geneIndexs, numOfClusters = 3, normalization = "normalizationZScore", discretization = "discretizationZScore",plotBase = ".")

# 	discretized_z = apply(sourceData, MARGIN = 1, discretizationZScore, numOfStates)
# 	oneDJuryScores_z = oneDJuryScore(discretized_z, geneIndexs, numOfStates = numOfStates)
# 	distances = dist(oneDJuryScores_z)
# 	fit_oneD_z = hclust(distances)
# 	groups_oneD_z = cutree(fit_oneD_z, k = numOfClusters)
# 	results_oneD_z = RRand(groups_oneD_z, trueLabel)$adjRand


# 	discretized_i = apply(sourceData, MARGIN = 1, discretizationFloor, numOfStates)
# 	oneDJuryScores_i = oneDJuryScore(discretized_i, geneIndexs, numOfStates = numOfStates)
# 	distances = dist(oneDJuryScores_i)
# 	fit_oneD_i = hclust(distances)
# 	groups_oneD_i = cutree(fit_oneD_i, k = numOfClusters)
# 	results_oneD_i = RRand(groups_oneD_i, trueLabel)$adjRand

# 	discretized_q = apply(sourceData, MARGIN = 1, discretizationQuantile, numOfStates)
# 	oneDJuryScores_q = oneDJuryScore(discretized_q, geneIndexs, numOfStates = numOfStates)
# 	distances = dist(oneDJuryScores_q)
# 	fit_oneD_q = hclust(distances)
# 	groups_oneD_q = cutree(fit_oneD_q, k = numOfClusters)
# 	results_oneD_q = RRand(groups_oneD_q, trueLabel)$adjRand

# 	against_trueLabel = data.frame(original = results_sourceData, quantile = results_Quantile, interval = results_Floor, zScore = results_zScore,
# 						oneD_feq = results_oneD_q, oneD_int = results_oneD_i, oneD_z = results_oneD_z)

# 	against_noOneD = data.frame(feq = RRand(groups_Quantile, groups_oneD_q)$adjRand,
# 								interval = RRand(groups_Floor, groups_oneD_i)$adjRand,
# 								zScore = RRand(groups_zScore, groups_oneD_z)$adjRand)
# 	result = list(vsTrueLabel = against_trueLabel, against_noOneD = against_noOneD)
# 	save(file = file.path("./processedData/oneDJury/", paste(c(fileName, numOfStates), collapse = "_")), result)
# 	result
# }