#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)

#read file
sData = read.csv("./pam50_label_orginal.xls", sep = "\t", header=TRUE);
#process header information(for labeling)
hInfo = colnames(sData);
hInfo = unlist(lapply(hInfo, processHeader));

# #test the performance of different scaling, distance and clustering functions 
# sourceData = sData[,hInfo > 0];
# sourceData = as.matrix(sourceData);
# benchmark = hInfo[hInfo > 0];
# savePath = "./plots/pam50/continuous_with_NA"
# plotFunctions(sourceData, savePath);
# #test the performance of different scaling, distance and clustering functions, exclude the NA group
# sourceData = sData[,hInfo > 0 & hInfo < 6];
# sourceData = as.matrix(sourceData);
# benchmark = hInfo[hInfo > 0 & hInfo < 6];
# savePath = "./plots/pam50/continuous_WO_NA"
# plotFunctions(sourceData, savePath);

# sourceData = sData[,hInfo > 0];
# sourceData = as.matrix(sourceData);
# benchmark = hInfo[hInfo > 0];
# savePath = "./plots/pam50/discrete_with_NA"
# plotStates(sourceData, savePath, 6, benchmark);

sourceData = sData[,hInfo > 0 & hInfo < 5];
sourceData = as.matrix(sourceData);
benchmark = hInfo[hInfo > 0 & hInfo < 5];
basePath = "./plots/pam50"
states  = 3  : 10
normalization = c("normalizationZScore","normalizationLinear")
discretization = c("discretizationZScore", "discretizationFloor")
for(i in 1 : length(discretization)) {
	titleName = discretization[i];
	bLine = baseLine(sourceData, numOfClusters = 4, normalization = match.fun(normalization[i]));
	similarity = pairWiseComparision(sourceData, numOfClusters = 4,
		benchMark = benchmark,states = states, discretization = match.fun(discretization[i]))
	plotEvaluations(similarity, states = states, baseLine = bLine,
	 	titleName = titleName, savePath = basePath, lx = 7.2, ly = 1)
}