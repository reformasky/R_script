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


evaluateOneDJury = function(sData, genesInGroups, normalization,  plotBase, discretization, numOfClusters,
	dendro = 3,  numOfStates = 3 : 5) {
	titleName = paste(c(normalization, "noDiscretization","dengro"), collapse = "_")
	trueLabeling = as.numeric(colnames(sData));
	sourceData = apply(sData, MARGIN = 1, match.fun(normalization));
	distance = dist(sourceData);
	fit = hclust(distance)
	nonDiscretizedGroups = cutree(fit, k = numOfClusters)
	# plotDendrogram(fit, savePath = plotBase, 
	# 	titleName = titleName);

	euclideanVsTrueLable = replicate(length(numOfStates), 0);
	euclieanVsNoDiscretize = replicate(length(numOfStates), 0);


	for(i in 1 : length(numOfStates)) {
		state = numOfStates[i]
		sourceData = oneDJuryScore(sData = sData, genesInGroups = genesInGroups, numOfStates = state, discretization = discretization);
		distance = dist(sourceData);
		fit = hclust(distance)
		if(state == dendro){
			# plotDendrogram(fit, savePath = plotBase, 
			# 	title = paste(c(discretization, "numOfStates", state), collapse = "_"))
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