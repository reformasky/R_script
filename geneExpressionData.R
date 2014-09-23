setwd("D:/thesis")
source("./R_script/commons.R");
fileNamesList = c("pam50", "pam50_p53", "miller_er", "miller_P53DLDA","miller_elston");
sampleTypesList = list(c("Basal","Her2","LumA","LumB"),
					c("wild", "mutated"),
					c("ER"),
					c("X0","X1"),
					c("G1", "G2", "G3"));
numOfClustersList = c(4, 2, 2, 2,3);
statesList = list(3 : 10, 3 :5, 3: 5, 3 : 5, 3 : 5)

for( i in 1 : length(fileNamesList)) {
	geneExpressionData(sourceFile = fileNamesList[i], sampleTypes = unlist(sampleTypesList[i]), numOfClusters = numOfClustersList[i], states = unlist(statesList[i]))
}


