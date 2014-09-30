source("./R_script/util.R")

#return a list of length numOfGroups, 
generateGroups = function(numOfGenes, numOfGroups) {
	genesInGroups = sample(1 : numOfGenes);
	genesPerGroup = numOfGenes / numOfGroups;
	result = list();
	for(g in 1 : numOfGroups) {
		result[g] = list(genesInGroups[((g - 1) * genesPerGroup + 1) : (g * genesPerGroup)])
	}
	result;
}


# sData: matrix of expression without normalization of discretization; genes x samples;
# genesInGroups(very important): a list of lists of genes, each inner list is the geneIds 
	# for that particular groups of genes; need to use unlist to retrive an array of genes
# returns the numOfSample x numOfGroups matrix, each row is a 1DJury score vector.
oneDJuryScore = function(sData, genesInGroups, numOfStates = 3, discretization = "discretizationZScore") {
	# assuming mtx is a matrix with proper discretization; returns nStates x numOfTotalGenes;
	# mtx is transposed because of discretization! numOfSamples x numOfGenes
	oneDJuryCache = function(mtx, numOfStates) {
		result = matrix(0, ncol = dim(mtx)[2], nrow = numOfStates);
		rownames(result) = 1 : numOfStates - 1;
		colnames(result) = colnames(mtx);
		for(r in 1: dim(result)[1]) {
			for(c in 1 : dim(result)[2]) {
				result[r, c] = sum(mtx[,c] == r );
			}
		}
		result;
	}

	numOfGroups = length(genesInGroups);

	#this is transposed !
	sourceData = apply(sData, MARGIN = 1, match.fun(discretization), numOfStates);
	#add 1 to avoid the indexing problem
	sourceData = sourceData + 1;
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
		result[,group] = t(rowSums(cachedSamples[, selectedGenes]))
	}
	
	result;
}


evaluateOneDJury = function(sData, genesInGroups, normalization,  discretization, numOfClusters, dendro = 3, plotBase ) {
	titleName = paste(c(normalization, "noDiscretization","dengro"), collapse = "_")
	trueLabeling = as.numeric(colnames(sData));
	sourceData = apply(sData, MARGIN = 1, match.fun(normalization));
	distance = dist(sourceData);
	fit = hclust(distance)
	nonDiscretizedGroups = cutree(fit, k = numOfClusters)
	plotDendrogram(fit, savePath = plotBase, 
		titleName = titleName);
	randAgainstTrueLabel = replicate(length(numOfStates), 0);
	randAgainstNoDiscretization = replicate(length(numOfStates), 0);
	for(i in 1 : length(numOfStates)) {
		state = numOfStates[i]
		sourceData = oneDJuryScore(sData = sData, genesInGroups = genesInGroups, numOfStates = state, discretization = discretization);
		distance = dist(sourceData);
		fit = hclust(distance)
		if(state == dendro){
			plotDendrogram(fit, savePath = plotBase, 
				title = paste(c(discretization, "numOfStates", state), collapse = "_"))
		}
		groups = cutree(fit, k = numOfClusters);
		randAgainstTrueLabel[i] = RRand(trueLabeling, groups)$adjRand
		randAgainstNoDiscretization[i] = RRand(nonDiscretizedGroups, groups)$adjRand;		
	}
	result = data.frame( vsTrueLabel = randAgainstTrueLabel, vsNoDiscretization = randAgainstNoDiscretization)
	rownames(result) = numOfStates;
	result = format(result, digits = 3);
	result = list(discretized = result, noDiscretizeVsTrueLabel = RRand(trueLabeling, nonDiscretizedGroups)$adjRand)
}