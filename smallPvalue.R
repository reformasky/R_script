#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

sData = read.csv("./TVGA_BreastCancer_fold3_p0.00001.xls", sep = "\t", header = TRUE)
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader))
colnames(sData) = hInfo;
#sort basing on p value 
sData = sData[order(sData$"0"),]

numOfGenes = 1 : 10 * 20;

basePath = "./plots/pvalues"

for(nGenes in numOfGenes) {
	titleName = paste(c("Number of genes_", nGenes), collapse = "");
	sourceData = sData[1: nGenes, hInfo > 0 & hInfo < 5];
	benchmark = hInfo[hInfo > 0 & hInfo < 5];
	similarity = pairWiseComparision(sourceData, savePath = basePath, 4, benchMark = benchmark,states = 3 : 5, 
	titleName = titleName, lx = 4.2, ly = 0.9)
}