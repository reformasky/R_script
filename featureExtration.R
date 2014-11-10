setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")

plotBase = "plots/featureExtraction"
saveBase = "processedData/featureExtraction"
discretizations = c("discretizationQuantile", "discretizationFloor", "discretizationZScore")
numOfStates = 3
numOfClusters = 4

objectBase = "processedData/oneDJury/TCGA_4Cancers_savedObjects"
load(file.path(objectBase, "combinedMatrix"))


zeroMeanOneSd = function(mtx) {	
	means = apply(mtx, MARGIN =2, mean)
	subtract = function(vec, means) {
		vec - means;
	}
	mtx = t(apply(mtx, MARGIN =1, subtract, means))

	suqardSum  = apply(mtx, MARGIN =2,  function(vec){sum(vec **2)})

	divide = function(vec, suqardSum) {
		vec / sqrt(suqardSum)
	}

	mtx = t(apply(mtx, MARGIN = 1, divide, suqardSum))
	mtx = mtx * sqrt(dim(mtx)[1] -1)
}

cNames = colnames(combinedMatrix)
rNames = rownames(combinedMatrix)
combinedMatrix = zeroMeanOneSd(combinedMatrix)
colnames(combinedMatrix) = cNames;
rownames(combinedMatrix) = rNames;


geneListFileName = "GPLists.csv"
geneLists = gpLists(geneListFileName)
combinedMatrix = findCombinedList(combinedMatrix, geneLists)
geneIndexs = findIndexesFromCombined(combinedMatrix, geneLists)

groups_linear = cutree(fit_linear, k = numOfClusters)
groups_zscore = cutree(fit_zscore, k = numOfClusters)
trueLabel = as.numeric(colnames(combinedMatrix))

normalized_linear = apply(combinedMatrix, MARGIN = 1, normalizationLinear)
distances = dist(normalized_linear)
fit = hclust(distances)
groups_linear = cutree(fit, k = numOfClusters)


normalized_zscore = apply(combinedMatrix, MARGIN = 1, normalizationZScore)
distances = dist(normalized_zscore)
fit = hclust(distances)
groups_zscore = cutree(fit, k = numOfClusters)



tops = c(5, 10, 20, length(geneIndexs))
result_vsNonDiscretizaion = replicate(length(tops), 0)
result_trueLabel = replicate(length(tops), 0)

result = c()

i = 1

for(i in 1 : length(discretizations)) {


	discretization = discretizations[i]

	sourceData = apply(combinedMatrix, MARGIN = 1, discretization, numOfStates = numOfStates)
	sourceData = oneDJuryScore(sData = sourceData, genesInGroups = geneIndexs, numOfStates = numOfStates)
	rownames(sourceData) = as.numeric(colnames(combinedMatrix))

	for(j in 1 : length(tops)) {
		selectedData = selectForTop(sData = sourceData, numOfGpsKept = tops[j])
		distance = dist(selectedData)
		fit = hclust(distance)
		groups = cutree(fit, k = numOfClusters)
		result_trueLabel[j] = RRand(groups, trueLabel)$adjRand
		if(i <= 2) {
			result_vsNonDiscretizaion[j] = RRand(groups, groups_linear)$adjRand
		}
		else {
			result_vsNonDiscretizaion[j] = RRand(groups, groups_zscore)$adjRand
		}

		titleName = paste(c(discretization, tops[j]), collapse = "_")
		plotDendrogram(fit = fit, savePath = plotBase, titleName = titleName)
	}

	if(i <=2) {
		nonDiscretizationVsTrueLabel = RRand(groups_linear, trueLabel)$adjRand
	}	else {
		nonDiscretizationVsTrueLabel = RRand(groups_zscore, trueLabel)$adjRand
	}		
	tempResult = list(data.frame(vsNonDiscretization = result_vsNonDiscretizaion, vsTrueLabel = result_trueLabel), nonDiscretizationVsTrueLabel)

	result = c(result, tempResult)
}

write.csv(file = file.path(saveBase, "featureExtraction.csv"), result)

tiff(file.path(plotBase, "featureExtraction.tiff"), units = "in", height = 6, width = 8, res = 300)
colors = c("red","green", "blue","pink")
par(mfcol  = c(2, length(discretizations)), mar = c(2,5,1,1))
for(i in 1 : length(discretizations)) {
	vsNonDiscretization = cbind(unlist(result[i *2 -1] )[1: 4])
	vsTrueLabel = cbind(unlist(result[i * 2 -1])[5: 8])
	referenceLine = unlist(result[i * 2])
	if( i == 1) {
		bp_non = barplot(vsNonDiscretization,beside = T, col = colors, xaxt = "n",  ylim = c(-0.05,1), ylab = "adjusted Rand Index")
		bp_true = barplot(vsTrueLabel,beside = T, col = colors, xaxt = "n",  ylim = c(-0.05,1), ylab = "adjusted Rand Index")
	} else {
		bp_non = barplot(vsNonDiscretization,beside = T, col = colors, xaxt = "n",  ylim = c(-0.05,1), yaxt = "n")
		bp_true = barplot(vsTrueLabel,beside = T, col = colors, xaxt = "n",  ylim = c(-0.05,1), yaxt = "n")
	}
	abline(referenceLine, 0, col = "black", lty = 3, lwd = 2)

}

dev.off()