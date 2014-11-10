setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")
saveBase = "./processedData/oneDJury/50Each"
plotBase = "./plots/oneDJury/TCGA_4Cancers/50Each"
load(file="combinedMatrix")

####
#### to select 50 sample from each group, remove comments
####

# colname = as.numeric(colnames(combinedMatrix))
# sampleSize = c(66, 392 + 66, 66+392+158, 66+392+158+ 246)
# cummulated = c(0,66,66 + 392, 66+392+158) + 1
# selected = c()
# for(i in 1 : 4) {
# 	selected = c(selected, sample(cummulated[i]:sampleSize[i], size =50))
# }
# combinedMatrix = combinedMatrix[,selected]
# colnames(combinedMatrix) = colname[selected]

numOfClusters = 4

normalizations = c("normalizationLinear", "normalizationLinear", "normalizationZScore")
discretizations = c("discretizationQuantile", "discretizationFloor", "discretizationZScore")
numOfStates = 3 : 5

geneListFileName = "GPLists.csv"
geneLists = gpLists(geneListFileName)
combinedMatrix = findCombinedList(combinedMatrix, geneLists)
geneIndexs = findIndexesFromCombined(combinedMatrix, geneLists)

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
# fileName = file.path(saveBase, paste(c("tcga_combined", "csv"), collapse = "."))
# write.csv(savedCSV, file = fileName)
# plotOneDJuryBar("tcga_combined")