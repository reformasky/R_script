#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./TheraPrint.csv", sep = ",", header = TRUE)
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader))
colnames(sData) = hInfo;

basePath = "./plots/theraPrint"

sourceData = sData[, hInfo > 0 & hInfo < 5];
benchmark = hInfo[hInfo > 0 & hInfo < 5];
numOfClusters = 3 : 10
for(nc in numOfClusters){
	titleName = paste(c("TheraPrint_numOfCLusters_", nc), collapse = "")

	similarity = pairWiseComparision(sourceData, savePath = basePath, numOfClusters = nc, benchMark = benchmark,states = 3 : 5, 
	titleName = titleName, lx = 4.2, ly = 0.9)
}