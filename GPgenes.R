setwd("D:/thesis")
gp = function(fileName) {
  a = file(fileName, open = "r")
  b = readLines(a)
  results = list()
  for( i in 1 : length(b)) {
    tempList = unlist(strsplit(b[i], ","))
    result = tempList[ tempList != "" ]
    results = c(results, list(result))
  }
  results
}

findIndexes = function(vec, lists) {
  results = list();
  for( i in 1 : length(lists) ) {
    lhs = unlist(lists[i])
    temp = match(lhs, vec);
    result = temp[!is.na(temp)]
    results = c(results, list(result))
  }
  results
}

