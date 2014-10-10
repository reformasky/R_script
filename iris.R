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
numOfClusters = 3

plotBarGraph(sourceData, 3, "iris", savePath = basePath)