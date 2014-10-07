setwd("D:/thesis")
source("./R_script/util.R")
savePath = "././processedData/discretizations"
plotBase = "./plots/extra_pam50"
set.seed(1024)

clusterings = c( "complete", "single", "ward")
normalizations = c("normalizationZScore","normalizationLinear", "normalizationLinear")
discretizations = c("discretizationZScore", "discretizationFloor", "discretizationQuantile")

sData = read.csv("pam50", sep = "\t", header=TRUE)
sampleTypes = c("Basal","Her2","LumA","LumB")
hInfo = colnames(sData)
hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));
sData = sData[,hInfo > 0 & hInfo <= length(sampleTypes)];
sData = as.matrix(sData);
benchmark = hInfo[hInfo > 0 & hInfo  <= length(sampleTypes)];
dendro = 3
colnames(sData) = as.numeric(benchmark)

numOfStates = 3
numOfClusters = 4



plotExtra = function(result, savePath, titleName,  baseline = 0, noDiscretizedVsTrueLabels = c()) {
	resultLength = dim(result)[1]
	color = "black"
	xLab = "Clustering Methods"
	yLab = "adjusted Rand index"
	fileName = file.path(savePath, paste(c(titleName, ".tiff"), collapse = ""))
	tiff(fileName,units="in", width=6, height=6, res=300)
	plot(1 : dim(result)[2], result[1,], xlab = xLab, ylab = yLab, xaxt = "n",
				main = titleName,  ylim = c(-0.05, 1), pch = 1, cex = 1.5, color = color)
	for(i in 2 : resultLength) {
		if( i > (resultLength /2))
			color = "red"
		par(new = T)
		plot(1 : dim(result)[2],result[i,], xlab = xLab, ylab = yLab, ylim = c( - 0.05, 1), 
		xaxt = "n", yaxt="n", pch = i, cex = 1.5, col = color);
	}
	axis(1, at = 1 : dim(result)[2], labels = colnames(result))
	abline(baseline, 0, col = "blue", lty = 0)
	for(i in 1 : length(noDiscretizedVsTrueLabels)) {
		abline(noDiscretizedVsTrueLabels[i], 0, col = "red", lty = i + 1)
	}
	dev.off()

}

#reference line is always complete, euclidean distance
for(i in 1 : length(normalizations)) {

	normalization = normalizations[i]
	discretization = discretizations[i]

	sourceData = apply(sData, MARGIN = 1, normalization)
	
	distance = dist(sourceData)
	fit = hclust(distance)
	plotDendrogram(fit, savePath = plotBase, titleName = paste(c(normalization, "noDiscetization"), collapse = "_"))
	noDiscretizedGroup = cutree(fit, k = numOfClusters)

	shuffleData = apply(t(sourceData), MARGIN = 1, sample)
	distances = dist(shuffleData)
	fit = hclust(distances)
	
	shuffledNoDiscretizedGroup = cutree(fit, k = numOfClusters)

	baseline = RRand(noDiscretizedGroup, shuffledNoDiscretizedGroup)$adjRand
	
	
	noDiscretizedVsTrueLabels = replicate(length(clusterings), 0)
	for(j in 1 : length(clusterings)) {
		clustering  = clusterings[j]

		sourceData = apply(sData, MARGIN = 1, normalization)
	
		distance = dist(sourceData)
		fit = hclust(distance, method = clustering)
		plotDendrogram(fit, savePath = plotBase, titleName = paste(c(normalization, "noDiscetization", clustering), collapse = "_"))
		noDiscretizedGroup = cutree(fit, k = numOfClusters)

		noDiscretizedVsTrueLabels[j] = RRand(noDiscretizedGroup, benchmark)$adjRand

		sourceData = apply(sData, MARGIN = 1, discretization)
		distance = dist(sourceData)
		fit = hclust(distance, method= clustering)
		plotDendrogram(fit, savePath = plotBase, 
			titleName = paste(c(normalization, clustering,"discretization","euclidean" ), collapse = "_"));
		groups = cutree(fit, k = numOfClusters)

		euclideanVsNoDiscretize = RRand(groups, noDiscretizedGroup)$adjRand
		euclideanVsTrueLabel = RRand(groups, benchmark)$adjRand

		distance = hammingMatrix(sourceData)
		fit = hclust(distance, method = clustering)
		plotDendrogram(fit, savePath = plotBase, 
			titleName = paste(c(normalization, clustering,"discretization", "hamming" ), collapse = "_"));
		groups = cutree(fit, k = numOfClusters)
		hammingVsNoDiscretize = RRand(groups, noDiscretizedGroup)$adjRand
		hammingVsTrueLabel = RRand(groups, benchmark)$adjRand

		distance = cosineMatrix(sourceData)
		fit = hclust(distance, method = clustering)
		plotDendrogram(fit, savePath = plotBase, 
			titleName = paste(c(normalization, clustering,"discretization", "cosine" ), collapse = "_"));
		groups = cutree(fit, k = numOfClusters)
		cosineVsNoDiscretize = RRand(groups, noDiscretizedGroup)$adjRand
		cosineVsTrueLabel = RRand(groups, benchmark)$adjRand
		result = c(euclideanVsNoDiscretize,hammingVsNoDiscretize, cosineVsNoDiscretize,
						euclideanVsTrueLabel, hammingVsTrueLabel, cosineVsTrueLabel)
		if(j == 1) {
			results = c(result)
		}
		else{
			results = c(results, result)
		}
		print(result)
	}
	results  = matrix(results, nrow = 6, byrow = F)
	colnames(results) = clusterings
	rownames(results) = c("euclideanVsNoDiscretize","hammingVsNoDiscretize", "cosineVsNoDiscretize",
						"euclideanVsTrueLabel", "hammingVsTrueLabel", "cosineVsTrueLabel")
	plotExtra(results, plotBase, titleName = discretization, baseline, noDiscretizedVsTrueLabels)
	write.csv(results, file.path(savePath, paste(c("extra_pam50", discretization), collapse = "")))
	print("\n******************************\n")
	print(i)
}
