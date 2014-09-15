library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./iris.txt", sep = ",", header = FALSE)

basePath = "./plots/iris"
states = 3 : 5
sourceData = sData[, 1 : 4];
sourceData = t(sourceData)
numOfClusters = 2 : 5
normalization = c("normalizationZScore","normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor")
benchMark = c(replicate(50, 1), replicate(50, 2), replicate(50,3))
for(i in 1 : length(discretization)) {
	for(nc in numOfClusters){
		titleName = paste(c(discretization[i],"iris_numOfCLusters", nc), collapse = "_")
		bLine = baseLine(sourceData, nc, normalization = match.fun(normalization[i]))
		similarity = pairWiseComparision(sourceData, numOfClusters = nc, 
			states = states, discretization = match.fun(discretization[i]), benchMark = benchMark);
		plotEvaluations(similarity, states, baseLine = bLine, 
			titleName = titleName, savePath = basePath, lx = 4.2, ly = 0.92)
	}
}