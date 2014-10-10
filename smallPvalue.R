#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./smallPValues.xls", sep = "\t", header = TRUE)
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader,sampleTypes = c("Basal","Her2","LumA","LumB")))
colnames(sData) = hInfo;
#sort basing on p value 
sData = sData[order(sData$"0"),]

numOfGenes = 1 : 10 * 20;

basePath = "./plots/smallPValues"
states = 3 : 5

for(nGenes in numOfGenes) {
	titleName = paste(c("numOfGenes",nGenes), collapse="_");
	sourceData = sData[1: nGenes, hInfo > 0 & hInfo < 5];
	colnames(sourceData) = hInfo[hInfo > 0 & hInfo < 5];
	plotBarGraph(sourceData, 4, titleName, savePath = basePath)
	print(nGenes)
}