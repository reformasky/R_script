setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./iris.txt", sep = ",", header = FALSE)
basePath = "./plots/iris"
savePath = "./processedData/discretizations"
states = 3 : 5
sourceData = sData[, 1 : 4];
sourceData = t(sourceData)
colnames(sourceData) = c(replicate(50, 1), replicate(50,2), replicate(50, 3))
normalization = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")
numOfClusters = 3
benchMark = as.numeric(colnames(sourceData))

for(i in 1 : length(discretization)) {
	titleName = paste(c(discretization[i],"iris"), collapse = "_")
	bLine = baseLine(sourceData, numOfClusters, normalization = match.fun(normalization[i]))
	results = pairWiseComparision(sourceData, numOfClusters = numOfClusters, states = states, 
		normalization = normalization[i], discretization = match.fun(discretization[i]), benchMark = benchMark, 
		dendro = 3, savePath = basePath,titleName = paste(c(titleName,"dendrograph"), collapse = "_"));
	plotEvaluations(results, states, baseLine = bLine, 
		titleName = titleName, savePath = basePath, lx = 4.2, ly = 0.92)
	similarity = data.frame(result$discretized)
		noDiscretizeVsTrueLabel = result$noDiscretizeVsTrueLabel
	colnames(similarity) = paste(discretization[i], colnames(similarity), sep = " ")
	if(i == 1) {
			savedData = data.frame(similarity);
			noDiscretization = data.frame(noDiscretizeVsTrueLabel)
		}
		else {
			savedData = data.frame(savedData, similarity);
			noDiscretization = data.frame(noDiscretization, noDiscretizeVsTrueLabel )
		}
}

write.csv(savedData, file = paste(c(savePath, "csv"), collapse = "."))
write.csv(noDiscretization, file = paste(c(savePath, "_noD", ".csv"), collapse = ""))