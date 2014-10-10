setwd("D:/thesis")
source("./R_script/oneDJuryUtil.R")

set.seed(1024)

numOfStates1  = 2

data1 = data.frame(a = c(0,0,1,0,1,0,0,0,1,1),
					b = c(0,0,0,0,1,0,0,0,1,1),
					c = c(1,0,0,1,0,1,1,1,0,0),
					d = c(0,1,1,1,1,1,0,1,1,0))
genesInGroups1 = c(list(1 : 10), list(1 : 9))
result1 = oneDJuryScore(data1, genesInGroups, 2)


numOfStates2 = 3
data2 = data.frame(a = sample(0:2, size = 5, replace = T),
				   b = sample(0:2, size = 5, replace = T),
				   c = sample(0:2, size = 5, replace = T),
				   d = sample(0:2, size = 5, replace = T), 
				   e = sample(0:2, size = 5, replace = T))

genesInGroups2 = c(list(c(2,5)), list(c(1,4)), list(1 : 5))
result2 = oneDJuryScore(data2, genesInGroups2, numOfStates = numOfStates2)

numOfStates3 = 3
data3 = data.frame(a = sample(0:2, size = 5, replace = T),
				   b = sample(0:2, size = 5, replace = T),
				   c = sample(0:2, size = 5, replace = T),
				   d = sample(0:2, size = 5, replace = T), 
				   e = sample(0:2, size = 5, replace = T),
				   f = sample(0:2, size = 5, replace = T))

genesInGroups3 = c(list(c(2,5)), list(c(1,4)), list(1 : 5))
result3 = oneDJuryScore(data3, genesInGroups3, numOfStates = numOfStates3)
