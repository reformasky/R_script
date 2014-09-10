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

sourceData = sData[,hInfo > 0];
sourceData = as.matrix(sourceData);
benchmark = hInfo[hInfo > 0];
savePath = "./plots/pam50/discrete_with_NA"
plotStates(sourceData, savePath, 6, benchmark);

sourceData = sData[,hInfo > 0 & hInfo < 6];
sourceData = as.matrix(sourceData);
benchmark = hInfo[hInfo > 0 & hInfo < 6];
savePath = "./plots/pam50/discrete_WO_NA"
plotStates(sourceData, savePath, 6, benchmark);