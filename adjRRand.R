
setwd("D:/thesis")
source("./R_script/util.R")
set.seed(1024)
numOfClusters = 4
numOfSamples = 200
numOfFeatures = 100

plotBase = "plots/forRandIndex"
trueLabel = sample(1 : 4, size = numOfSamples, replace = T)
centers = matrix(0, ncol = numOfClusters, nrow = numOfFeatures)
for(i in 1 : dim(centers)[1]) {
	centers[i,] = sample(1 : 4, size = numOfClusters, replace = T)
}

generateData = function(stdev) {
	result = matrix(0, ncol = numOfSamples, nrow = numOfSamples)
	for(i in 1 : numOfSamples) {
		result[,i] = rnorm(numOfFeatures, mean = centers[,trueLabel[i]], sd = stdev);
	}
	result;
}

stdevs = c(1.5, 2.1, 2.5,  2.65, 3, 3.27, 3.45, 3.51,  4.3,  7)

adjRand = 1 : length(stdevs)
rand = 1 : length(stdevs)
for(i in 1 : length(stdevs)) {
	set.seed(1024)
	stdev = stdevs[i]
	simulatedData = generateData(stdev);
	colnames(simulatedData) = trueLabel;
	mtx = t(simulatedData)
	fit = hclust(dist(mtx))
	groups = cutree(fit, k = numOfClusters)
	randIndexs = RRand(groups, trueLabel)
	adjRand[i] = randIndexs$adjRand
	rand[i] = randIndexs$Rand
	titleName = paste(c(i, "adjRand", randIndexs$adjRand, "Rand", randIndexs$Rand), collapse = "_")
	plotDendrogram(fit, titleName = titleName, savePath = plotBase)

}

write.csv(data.frame(adjRand, rand), file = file.path(plotBase, "numbers.csv"))





























# trueLabel = c(replicate(100, 1), replicate(100,2), replicate(100, 3), replicate(100,4))
# testLabel = replicate(400, 0)
# reservedNumber = c(1, 33, 44,52, 63, 73, 77,84, 90,96, 100)
# resultAdj = replicate(length(reservedNumber),0)
# resultRand = replicate(length(reservedNumber),0)

# for( index in 1 : length(reservedNumber)) {
# 	i = reservedNumber[index]
# 	reserved = c( 1: i, 1 : i + 100, 1 : i + 200, 1 : i + 300 )
# 	testLabel[reserved] = trueLabel[reserved]
# 	testLabel[-reserved] = sample(1 : 4, replace =T, size = length(trueLabel) - length(reserved))
# 	df = data.frame(seq = 1 : length(testLabel), cluster = testLabel);

# 	resultAdj[index] = RRand(trueLabel, testLabel)$adjRand
# 	resultRand[index] = RRand(trueLabel, testLabel)$Rand

# 	p2<-ggplot(df,aes(seq,y=1,fill=factor(cluster)))+ geom_tile()+
# 	  scale_y_continuous(expand=c(0,0))+
# 	  theme(axis.title=element_blank(),
# 	        axis.ticks=element_blank(),
# 	        axis.text=element_blank(),
# 	        plot.background = element_blank(),
#    			panel.grid.major = element_blank(),
#    			panel.grid.minor = element_blank(),
#    			panel.border = element_blank(),
# 	        legend.position="none");
# 	tiff(paste(c(i, ".tiff"), collapse = ""))
# 	gp2<-ggplotGrob(p2)
# 	gp2$layout$clip[gp2$layout$name=="panel"] <- "off"
# 	grid.arrange(gp2, ncol=1); 
# 	dev.off()
#  }

#  result = data.frame(resultAdj, resultRand);
#  colnames(result) = c("adjRand","Rand");
#  rownames(result) = reservedNumber;
#  write.csv(result, file = "number.csv");