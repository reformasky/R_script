setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")
set.seed(1024)

normalizations = c("normalizationLinear", "normalizationLinear", "normalizationZScore")
discretizations = c("discretizationQuantile", "discretizationFloor", "discretizationZScore")
numOfStates = 3 : 5

plotBase = "./plots/examineTCGA"
saveBase = "processedData/tcga_examine"
load("processedData/oneDJury/TCGA_4Cancers_savedObjects/combinedMatrix")

#zeroMeanOneSd
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

plotBase = "./plots/examineTCGA/zeroMeanOneSd"
load("processedData/oneDJury/TCGA_4Cancers_savedObjects/combinedMatrix")
rNames = rownames(combinedMatrix)
combinedMatrix = zeroMeanOneSd(combinedMatrix)
quart1 = apply(combinedMatrix, MARGIN = 2, quantile, 0.25)
means = apply(combinedMatrix, MARGIN = 2, mean)
quart3 = apply(combinedMatrix, MARGIN = 2, quantile, 0.75)
sds = apply(combinedMatrix, MARGIN = 2, sd)
tiff(file.path(plotBase,"combined_box_zeroMeanOneSd.tiff"), units="in", width=6, height=12, res=300)
par(mfcol  = c(4, 1))
boxplot(quart1 ~ colnames(combinedMatrix))
title("0.25 quantile")
boxplot(means ~ colnames(combinedMatrix))
title("mean")
boxplot(quart3 ~ colnames(combinedMatrix))
title("0.75 quantile")
boxplot(sds ~ colnames(combinedMatrix))
title("sds")
dev.off()
geneListFileName = "GPLists.csv"
geneLists = gpLists(geneListFileName)
combinedMatrix = findCombinedList(combinedMatrix, geneLists)
geneIndexs = findIndexesFromCombined(combinedMatrix, geneLists)
numOfClusters = 4
labels =as.numeric(colnames(combinedMatrix))

rownames(combinedMatrix) = rNames
for(j in 1 : length(discretizations)) {
		results = evaluateOneDJury(sData = combinedMatrix, genesInGroups = geneIndexs, numOfClusters = numOfClusters, normalization = normalizations[j],
		        discretization = discretizations[j], numOfStates = numOfStates, plotBase = plotBase,dendro = 3)
		if(j == 1) {
		savedCSV = data.frame(results)
		}
		else {
			savedCSV = data.frame(savedCSV,  results);
		}
		print(discretizations[j])
	}
fileName = file.path(saveBase, paste(c("tcga_combined_zeroMeanOneSd", "csv"), collapse = "."))
write.csv(savedCSV, file = fileName)
plotOneDJuryBar("tcga_combined_zeroMeanOneSd", dataBase = saveBase)