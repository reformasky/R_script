### librarys
library("cgdsr")
library(org.Hs.eg.db)
library(annotate)

setwd("D:/thesis")
set.seed(1024)

geneListFile = "GPLists.csv"
openedFile = readLines(file(geneListFile, open = "r"))
geneIds = c()
for(i in 1: length(openedFile)) {
	geneIds = union(geneIds, unlist(strsplit(openedFile[i], ",")))
}
geneIds = geneIds[!(geneIds == "")]
geneList = getSYMBOL(geneIds, data='org.Hs.eg')
#we need to discard 71 genes, so discard randomly
geneList = sample(geneList)

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
studiesIndex = c(28, 30,49,55)


#studies Id for chromphobe, rental clear cell, ovarian servious and prostate ad
studiesId = c("kich_tcga", "kirc_tcga_pub","ov_tcga","prad_tcga")
sampleSize = c(66, 392, 158, 246)
cummulated = c(0,66,66 + 392,66+392+158)

#download data
# for(i in 1 : length(studiesId)){

# 	print("**********")
# 	print(i)
# 	print("**********")
	
# 	# [1, 1] returns a unique ID used to identify the case list ID in subsequent interface calls
# 	mycaselist = getCaseLists(mycgds,studiesId[i])[1, 1]
# 	#[3,1] returns 
# 	mygeneticprofile = getGeneticProfiles(mycgds,studiesId[i])[3,1]

# 	for(j in 1 : 68){
# 		if(j == 1)
# 		downloadData = getProfileData(mycgds,geneList[(j-1) * 100 + 1 : 100],mygeneticprofile,mycaselist)
# 		else
# 		downloadData = cbind(downloadData, getProfileData(mycgds,geneList[(j-1) * 100 + 1 : 100],mygeneticprofile,mycaselist))

# 		print(j)
# 	}
# 	#remove duplicated data? "XK"
# 	downloadData = downloadData[, colnames(downloadData) != "XK"]
# 	downloadData = t(downloadData)

# 	#remove na, nan, inf
# 	downloadData = downloadData[complete.cases(downloadData * 0), ,drop = F]

# 	#saveDowloaded File
# 	save(file = studiesId[i], downloadData)
# 	colnames(downloadData) =  replicate(dim(downloadData)[2], i)
# 	downloadData = data.frame(id = rownames(downloadData), downloadData)

# 	if(i == 1) {
# 		combinedMatrix = downloadData
# 	}
# 	else {
# 		combinedMatrix = merge(combinedMatrix, downloadData, by = "id")
		
# 	}
# }



for(i in 1 : 4) {
	load(studiesId[i])
	colnames(downloadData) =  replicate(dim(downloadData)[2], i)
	downloadData = data.frame(id = rownames(downloadData), downloadData)

	if(i == 1) {
		combinedMatrix = downloadData
	}
	else {
		combinedMatrix = merge(combinedMatrix, downloadData, by = "id")
		
	}
}
geneHUGOs = combinedMatrix[,1]


matchedIndex = match(geneHUGOs, geneList)
naIndex = is.na(matchedIndex)
geneHUGOs = geneHUGOs[!naIndex]

matchedIndex = match(geneHUGOs, geneList)
combinedMatrix = combinedMatrix[!naIndex, ]
rownames(combinedMatrix) = names(geneList)[matchedIndex]
combinedMatrix = combinedMatrix[,2 :dim(combinedMatrix)[2]]
colnames(combinedMatrix) = c(replicate(sampleSize[1], 1), replicate(sampleSize[2], 2),replicate(sampleSize[3], 3),replicate(sampleSize[4], 4))
save(file = "combinedMatrix", combinedMatrix)
#select 50 genes from each group
selected = c()
for(i in 1 : 4) {
	selected = c(selected, sample((1 : sampleSize + cummulated[i]), size  = 50, replace = F))
}

selectedMatrix = combinedMatrix[, selected[c(1:50, 151:200)]]