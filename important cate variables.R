library(MASS)
library(cluster)

setwd("/Users/sundaycat/Documents/Math 285/project")

load("full_imputed_data.R")

# separate the categorical and continuous data
catIndx  <- unlist(lapply(1:78, function(index,o) class(o[,index]) =="factor", imp.train))
cateData <- imp.train[,catIndx]  # pure categorical data set
contData <- imp.train[,!catIndx] # pure continuous data set

# perform MCA on pure categorical data set
len <- length(cateData)
mca.out <- MCA(cateData, graph = FALSE, ncp = len - 1)

# decide the number of principle components based on the elbow
m.eig <- mca.out$eig$eigenvalue
m.eig <- ifelse(m.eig > 1/len, ((len/(len-1))*(m.eig-1/len))^2, 0)

plot(1:(len-1), m.eig[1:(len-1)], bty="n", pch=20)
lines(1:(len-1), m.eig[1:(len-1)])

# pick first 6 component to analysis
col <- c(1:6) 
dimContrib <- mca.out$var$contrib[,col]/100
#wtContrib <- data.frame(dimContrib %*% (m.eig/sum(m.eig))[col])
wtContrib <- data.frame(dimContrib %*% m.eig[col])

# remove the level index for each variables names
name <- rownames(wtContrib)
pos <- regexpr(pattern ="_[^_]*$", name)
name <-ifelse(pos == -1, name, substr(name,1,pos-1))

# create weighted contribution table
rawTable <- cbind(name, data.frame(dimContrib, wtContrib))
cateTable <- aggregate.data.frame(rawTable[,2:8], by=list(rawTable$name), FUN = sum)
colnames(cateTable) <- c("Name",paste("MCA", col), "contribution score")
cateTable <- cateTable[order(cateTable[,length(col)+2],decreasing = T), ]

# save contribition table to an excel file
cateTable[,2:8] <- round(cateTable[,2:8], 5)
write.csv(x = cateTable, file = "categorical contribution score.csv")

