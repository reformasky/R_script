library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./memmaPrint.csv", sep = "\t", header = TRUE)
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader))
colnames(sData) = hInfo;

basePath = "./plots/memmaPrint"
states = 3 : 5
sourceData = sData[, hInfo > 0 & hInfo < 5];
numOfClusters = 2 : 10
normalization = c("normalizationZScore","normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor")
for(i in 1 : length(discretization)) {
	for(nc in numOfClusters){
		titleName = paste(c(discretization[i],"memmaPrint_numOfCLusters", nc), collapse = "_")
		bLine = baseLine(sourceData, nc, normalization = match.fun(normalization[i]))
		similarity = pairWiseComparision(sourceData, numOfClusters = nc,states = states, 
			normalization = normalization[i],discretization = match.fun(discretization[i]));
		plotEvaluations(similarity, states, baseLine = bLine, 
			titleName = titleName, savePath = basePath, lx = 4.2, ly = 0.92)
	}
}