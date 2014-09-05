#library used to calculate adjusted Rand index
library(phyclust, quiet = TRUE)

setwd("D:/thesis")
source("./R_script/util.R")
#read file
sourceData = read.csv("./pam50_label_orginal.xls", sep = "\t", header=TRUE);

#process header information(for labeling)
headerInfo = colnames(sourceData);
headerInfo = unlist(lapply(headerInfo, processHeader));

sourceData = sourceData[,headerInfo > 0 ];
sourceData = as.matrix(sourceData);

benchmark = headerInfo[headerInfo > 0  ];

#parameters for hclust;
scalingMethodsString = c("zeroOneScaling", "zscoreScaling");
distanceMethods = c( "euclidean", "maximum", "manhattan",  "binary", "minkowski","canberra");
clusterMethods = c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid");

#zeroOneScaling
for(sMethod in scalingMethodsString){
	for(dMethod in distanceMethods) {
		for(cMethod in clusterMethods) {
			result = paste(c(sMethod, dMethod, cMethod), collapse = "_")
			print(paste(c(result, 
				evaluation(sourceData, 6, sMethod, dMethod, cMethod, benchmark)$Rand), collapse= ": "))
		}
	}
}