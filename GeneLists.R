setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")
saveBase = "./processedData/oneDJury"
fileNames = c("Collins_PAM50", "Collins_LUAD", "Collins_Colon")
numOfClustersList = c(4, 2, 2)
sampleTypesList = c(list(c("Basal","Her2","LumA","LumB")),
					list(c("A","B")),
					list(c("A", "B")))

normalizations = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretizations = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfStates = 3 : 5

for( i in 1 : length(fileNames))
{
	fileName = fileNames[i]
	numOfClusters = numOfClustersList[i]
	sampleTypes = unlist(sampleTypesList[i])
	dataFileName = file.path("data", paste(c(fileName,".txt"), collapse = ""))
	sData = read.csv(dataFileName, sep = "\t", header = T)
	hInfo = colnames(sData)
	hInfo = unlist( lapply(hInfo, processHeader, sampleTypes) )
	colnames(sData) = hInfo 

	sourceData = sData[, hInfo > 0]
	#only retains the first appearance of a gene id
	geneIds = unique(sData[,1])
	geneIndexs = match(geneIds, sData[,1])
	sourceData = sourceData[geneIndexs,]
	rownames(sourceData) = sData[geneIndexs,1]


	geneListFileName = "GPLists.csv"
	geneLists = gpLists(geneListFileName)
	sourceData = findCombinedList(sourceData, geneLists)
	geneIndexs = findIndexesFromCombined(sourceData, geneLists)

	trueLabel = hInfo[hInfo > 0]

	for(j in 1 : length(discretizations)) {
		results = evaluateOneDJury(sData = sourceData, genesInGroups = geneIndexs, numOfClusters = numOfClusters, normalizations[j], discretizations[j])
		if(j == 1) {
		savedCSV = data.frame(results)
		}
		else {
			savedCSV = data.frame(savedCSV,  results);
		}
		print(discretizations[j])
	}

	write.csv(savedCSV, file = file.path(saveBase, paste(c(fileName, "csv"), collapse = ".")))
	print(fileName)
}

for( i in 1 : length(fileNames))
{
	fileName = fileNames[i]
	plotOneDJuryBar(fileName)
}
	# distances = dist(t(sourceData))
	# fit_sourceData = hclust(distances)
	# groups_sourceData = cutree(fit_sourceData, k = numOfClusters)
	# results_sourceData = RRand(groups_sourceData, trueLabel)$adjRand

	# distances = dist(apply(sourceData, MARGIN = 1, discretizationZScore, 3))
	# fit_zScore = hclust(distances)
	# groups_zScore = cutree(fit_zScore, k = numOfClusters)
	# results_zScore = RRand(groups_zScore, trueLabel)$adjRand

	# distances = dist(apply(sourceData, MARGIN = 1, discretizationFloor, 3))
	# fit_Floor = hclust(distances)
	# groups_Floor = cutree(fit_Floor, k = numOfClusters)
	# results_Floor = RRand(groups_Floor, trueLabel)$adjRand

	# distances = dist(apply(sourceData, MARGIN = 1, discretizationQuantile, 3))
	# fit_Quantile = hclust(distances)
	# groups_Quantile = cutree(fit_Quantile, k = numOfClusters)
	# results_Quantile = RRand(groups_Quantile, trueLabel)$adjRand

	# #oneDJuryScores = oneDJuryScore(sourceData, geneIndexs)

	# #results = evaluateOneDJury(sourceData, geneIndexs, numOfClusters = 3, normalization = "normalizationZScore", discretization = "discretizationZScore",plotBase = ".")

	# discretized_z = apply(sourceData, MARGIN = 1, discretizationZScore, 3)
	# oneDJuryScores_z = oneDJuryScore(discretized_z, geneIndexs, numOfStates = 3)
	# distances = dist(oneDJuryScores_z)
	# fit_oneD_z = hclust(distances)
	# groups_oneD_z = cutree(fit_oneD_z, k = numOfClusters)
	# results_oneD_z = RRand(groups_oneD_z, trueLabel)$adjRand


	# discretized_i = apply(sourceData, MARGIN = 1, discretizationFloor, 3)
	# oneDJuryScores_i = oneDJuryScore(discretized_i, geneIndexs, numOfStates = 3)
	# distances = dist(oneDJuryScores_i)
	# fit_oneD_i = hclust(distances)
	# groups_oneD_i = cutree(fit_oneD_i, k = numOfClusters)
	# results_oneD_i = RRand(groups_oneD_i, trueLabel)$adjRand

	# discretized_q = apply(sourceData, MARGIN = 1, discretizationQuantile, 3)
	# oneDJuryScores_q = oneDJuryScore(discretized_q, geneIndexs, numOfStates = 3)
	# distances = dist(oneDJuryScores_q)
	# fit_oneD_q = hclust(distances)
	# groups_oneD_q = cutree(fit_oneD_q, k = numOfClusters)
	# results_oneD_q = RRand(groups_oneD_q, trueLabel)$adjRand

	# against_trueLabel = data.frame(original = results_sourceData, quantile = results_Quantile, interval = results_Floor, zScore = results_zScore,
	# 					oneD_feq = results_oneD_q, oneD_int = results_oneD_i, oneD_z = results_oneD_z)

	# against_noOneD = data.frame(feq = RRand(groups_Quantile, groups_oneD_q)$adjRand,
	# 							interval = RRand(groups_Floor, groups_oneD_i)$adjRand,
	# 							zScore = RRand(groups_zScore, groups_oneD_z)$adjRand)
	# result = list(vsTrueLabel = against_trueLabel, against_noOneD = against_noOneD)
	# save(file = file.path("./processedData/oneDJury/", fileName), result)
