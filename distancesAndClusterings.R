setwd("D:/thesis")

source("./R_script/util.R")
fileNamesList = c("pam50", "pam50_p53", "miller_er", "miller_P53DLDA","yipeng");
sampleTypesList = list(c("Basal","Her2","LumA","LumB"),
					c("wild", "mutated"),
					c("POS", "NEG"),
					c("X0","X1"),
					# c("G1", "G2", "G3"),
					c("quart1","quart4"));
numOfCLustersList =  c(4, 2, 2, 2, 4);
normalizationList = c("normalizationZScore", "normalizationLinear")


results = list()
for(i in 1 : length(fileNamesList)) {
	sampleTypes = unlist(sampleTypesList[i]);
	sourceFileName = fileNamesList[i];

	sData = read.csv(sourceFileName, sep = "\t", header=TRUE);
	hInfo = colnames(sData);
	hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));
	sourceData = sData[,hInfo > 0 & hInfo <= length(sampleTypes)];
	sourceData = as.matrix(sourceData);
	colnames(sourceData) = hInfo[hInfo > 0 & hInfo <= length(sampleTypes)];
	

	
	j = 1
	# for(j in 1 : length(normalizationList)) {
	titleName = paste(c(fileNamesList[i], normalizationList[j]), collapse = "_");
	result = evaluateDistanceAndClustering(sourceData, scalingMethod = normalizationList[j],numOfClusters = numOfCLustersList[i],titleName = titleName)
	# }	
	results = c(results, list(result))

}

tiff(file.path("./plots/distancesAndClusterings", paste(c("distances", ".tiff"), collapse = "")), units="in", width=8, height=5, res=300)
par(mfcol  = c(2, 1), mar = c(2,5,1,1))
colors = c("black", "red", "yellow", "green", "pink","blue")
fNames = c("Collins_PAM50", "Collins_p53", "Miller_ER", "Miller_p53DLDA", "miller_elston", "Wang_Stroma")

plotData1 = cbind(unlist(results[1])[c(1,3,5,7,9,11)])
plotData2 = cbind(unlist(results[1])[c(1,3,5,7,9,11) + 1])

for(i in 2 : length(fileNamesList)) {
	plotData1 = cbind(plotData1, unlist(results[i])[c(1,3,5,7,9,11)])
	plotData2 = cbind(plotData2, unlist(results[i])[c(1,3,5,7,9,11) + 1])
}

barplot(plotData1, beside = T, col = colors, main = "Euclidean", ylim = c(0, 1), ylab = "adjusted Rand Index")
barplot(plotData2, beside = T, col = colors, main = "Manhattan", ylim = c(0, 1), ylab = "adjusted Rand Index")
dev.off()




