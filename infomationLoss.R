setwd("D:/thesis")
fileNamesList = c("pam50", "pam50_p53", "miller_er", "miller_P53DLDA","yipeng");


for(i in 1 : length(fileNamesList)) {
	fileName = paste(c(fileNamesList[i], ".csv"), collapse = "")
	result = read.csv(file.path("processedData/discretizations", fileName), header  =T)
	e = c()
	h = c()
	for(j in 1: 3){
		e = c(e, (result[3,2 : 4 + (j -1) * 4] / result[3,j * 4 + 1]))
		h = c(h, (result[4,2 : 4 + (j -1) * 4] / result[4,j * 4 + 1]))
	}
	if(i == 1) {
		euclidean = data.frame(e)
		hamming = data.frame(h)
	}
	else{
		euclidean = data.frame(euclidean, e)
		hamming = data.frame(hamming, h)
	}
}
euclidean = as.numeric(euclidean)
hamming = as.numeric(hamming)
euclidean = matrix(euclidean, nrow = length(fileNamesList), byrow = T)
hamming = matrix(hamming, nrow = length(fileNamesList), byrow = T)
indexes = 1 :9
result = 1 : 9
for(i in 1 : 3){
	result[1 : 3 + (i -1)*3] = indexes[c(1, 4, 7) + i -1]
}

plotBase = "./plots/discretizations"
euclidean = euclidean[, result]
hamming = hamming[,result]	
tiff(file.path(plotBase,"combined.tiff"), units="in", width=8, height=6, res=300)
par(mfrow  = c(2, 1), mar = c(2,5,1,1))
euclideanMeans = apply(euclidean, MARGIN = 2, mean)
euclideanSDs = apply(euclidean, MARGIN = 2, sd)
hammingMeans = apply(hamming, MARGIN = 2, mean)
hammingSDs = apply(hamming, MARGIN = 2, sd)

euclideanMeans = cbind(euclideanMeans[1:3], euclideanMeans[4:6], euclideanMeans[7:9])
euclideanSDs = cbind(euclideanSDs[1:3], euclideanSDs[4:6], euclideanSDs[7:9])
hammingMeans = cbind(hammingMeans[1:3], hammingMeans[4:6], hammingMeans[7:9])
hammingSDs = cbind(hammingSDs[1:3], hammingSDs[4:6], hammingSDs[7:9])

ePlot = barplot(euclideanMeans, beside = T, col = c("red","green","blue"), ylim = c(0,2), xxaxt = "n")
segments(ePlot, euclideanMeans, ePlot, euclideanMeans + euclideanSDs, lwd = 3)
segments(ePlot - 0.1, euclideanMeans + euclideanSDs, ePlot + 0.1, euclideanMeans + euclideanSDs, lwd = 3)

ePlot = barplot(hammingMeans, beside = T, col = c("red","green","blue"), ylim = c(0,2), xxaxt = "n")
segments(ePlot, hammingMeans, ePlot, hammingMeans + hammingSDs, lwd = 3)
segments(ePlot - 0.1, hammingMeans + hammingSDs, ePlot + 0.1, hammingMeans + hammingSDs, lwd = 3)
dev.off()
