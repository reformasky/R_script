#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)
savePath = "./plots/pam50"
#read file
sourceData = read.csv("./pam50_label_orginal.xls", sep = "\t", header=TRUE);
#process header information(for labeling)
headerInfo = colnames(sourceData);
headerInfo = unlist(lapply(headerInfo, processHeader));

sourceData = sourceData[,headerInfo > 0 ];
sourceData = as.matrix(sourceData);
benchmark = headerInfo[headerInfo > 0  ];

#parameters for hclust;
scalingMethods = c("zeroOneScaling", "zscoreScaling");
distanceMethods = c( "euclidean", "maximum", "manhattan",  "binary", "minkowski","canberra");
clusterMethods = c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid");
clusterLabels = c("wa","si","co","av","mc",",me","ce")
labels = replicate( length(clusterMethods), "")
results = replicate( length(clusterMethods), 1)

for(s in 1 : length(scalingMethods)){
	for(d in 1 : length(distanceMethods)){
		cat(paste(replicate(length(clusterMethods), "*"), collapse = ""))
		print("")
		title = paste(c(scalingMethods[s], distanceMethods[d]), collapse = "_")
		counter = 1;
		for(cl in 1 : length(clusterMethods)) {
			label = clusterLabels[cl];
			result = evaluation(sourceData, 6,
			 	scalingMethods[s], distanceMethods[d],
			 	clusterMethods[cl], benchmark)$Rand;
			labels[counter] = label;
			results[counter] = result;
			counter = counter + 1;
			cat(sprintf("-"))
		}
		print("")
		fileName = paste(c(savePath, title), collapse = "/")
		fileName = paste(c(fileName, "tiff"), collapse = ".")
		tiff(fileName);
		plot(factor(labels), results, main= title, xlab= "agglomeration methods", ylab = "Rand index", ylim= c(0.3,0.7))
		dev.off();
	}
	
}
