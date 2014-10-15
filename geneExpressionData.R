setwd("D:/thesis")
source("./R_script/commons.R");
fileNamesList = c("pam50", "pam50_p53", "miller_er", "miller_P53DLDA","miller_elston", "yipeng");
sampleTypesList = list(c("Basal","Her2","LumA","LumB"),
					c("wild", "mutated"),
					c("POS", "NEG"),
					c("X0","X1"),
					c("G1", "G2", "G3"),
					c("quart1","quart4"));
numOfClustersList = c(4, 2, 2, 2,3, 2);
statesList = list(3 : 5, 3 :5, 3: 5, 3 : 5, 3 : 5, 3 :5)

for( i in 1 : length(fileNamesList)) {
	result =	geneExpressionData(sourceFile = fileNamesList[i], sampleTypes = unlist(sampleTypesList[i]), numOfClusters = numOfClustersList[i], states = unlist(statesList[i]))
	write.csv(file = file.path("./processedData/discretizations", paste(c(fileNamesList[i], ".csv"), collapse = "") ), result)
}


