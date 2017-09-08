############################################################
######## MATH 285 - SPRING 2017 DATA CHALLENGE PROJECT #####
########         HONGZHE LIU, LE PHAN      #################
############################################################

rm(list=ls())

#install.packages("mice")
#install.packages("gridExtra")
library(reshape2) #for boxplot
library(ggplot2) #for plot
library(gdata) #for load data
library(mice)
library(gridExtra)
library(reshape2)

# set the working directory
setwd("/Users/sundaycat/Desktop/Math 285/Project")

############################ FUNCTIONS START ###############################
# functions thant use to consolidate the data set that create by MICE pkg
consolDataSet <- function(mice.out, dataSet) {
  
  data <- dataSet
  data1 <- complete(mice.out, action = 1)
  data2 <- complete(mice.out, action = 2)
  data3 <- complete(mice.out, action = 3)
  data4 <- complete(mice.out, action = 4)
  data5 <- complete(mice.out, action = 5)
  #i=5;j=78
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
      if (!is.na(data[i, j]))
        next
      if (class(data[, j]) == "factor") {
        v1 <- as.integer(as.vector(data1[i, j]))
        v2 <- as.integer(as.vector(data2[i, j]))
        v3 <- as.integer(as.vector(data3[i, j]))
        v4 <- as.integer(as.vector(data4[i, j]))
        v5 <- as.integer(as.vector(data5[i, j]))
        v <- c(v1, v2, v3, v4, v5)
        
        maxFreq <- sort(table(v), decreasing = T)[1]
        data[i, j] <- as.integer(names(maxFreq))
        
      } else{
        v <- c(data1[i, j], data2[i, j], data3[i, j], data4[i, j], data5[i, j])
        data[i, j] <- mean(v)
      }
    }
  }
  
  return(data)
}

# Use to plot variables with regarding to each cluster
# output: cluster membership, data: noImputeData, varIndx: index of 2w,3m,12m variables
# plotType: 1 for validate variables, 2 for normal variables(apply to ordinal and continuous)
plotClust <- function(output, data, varIndx, graphName, plotType){
  
  names(output) <- NULL
  hasMean <- ifelse(plotType == 2, 1 ,0)
  
  gps <- length(unique(output))
  vars <- length(varIndx)
  
  # get the observation index in each cluster
  cluList <- list()
  for(i in 1:gps) cluList[[i]] <- which(output == i)
  
  # one row represent one cluster line
  cluLineMat <- matrix(0, nrow = gps+hasMean, ncol = vars)
  for(g in 1:gps){
    for(v in 1:vars){
      cluLineMat[g,v] <- mean(data[cluList[[g]],varIndx[v]], na.rm = T)
    }
  }
  
  
  if(plotType == 1){
    colnames(cluLineMat) <- c("2 weeks","3 months","12 months")
    rownames(cluLineMat) <- paste("Clust", 1:gps)
  }
  
  if(plotType == 2){
    cluLineMat[gps+hasMean, ] <- colMeans(as.matrix(data[,varIndx]), na.rm = T)
    colnames(cluLineMat) <- names(data)[varIndx]
    rownames(cluLineMat) <- c(paste("Clust", 1:gps),"Mean")
  }
  
  plotDf <- melt(cluLineMat)  #the function melt reshapes it from wide to long
  names(plotDf) <- c("Class", "Variable", graphName )
  
  return(plotDf)
}

# Ust to plot categorical Data
# data:noImputeData, cateVarName: name of variable, cluster: obs membership
#ordPlotData, "dlva0", k.cluster
plotCateVar <- function(noImputeData, cateVarName, cluster){
  
  varIndx <- which(names(noImputeData) == cateVarName)
  gps <- length(unique(cluster))
  numOflevel <- length(levels(noImputeData[,varIndx])) # how many levels
  levels <- as.integer(levels(noImputeData[,varIndx])) # possible values of levels
  # need to determined if there is NA
  isNa <- ifelse(anyNA(noImputeData[,varIndx]), 1 ,0)
  # create a matrix here
  pltMat <- matrix(0, nrow=numOflevel+isNa, ncol = gps)
  naCount <- vector(length=gps)
  for(g in 1 : gps){
    # some variables level start from 1, some starts from 0
    clustInx <- which(cluster == g)
    pltMat[1:numOflevel, g] <- sapply(levels, function(idx,o) length(which(o == idx)), noImputeData[clustInx,varIndx])
    naCount[g] <- length(which(is.na(noImputeData[clustInx,varIndx])))
  }
  if(isNa == 1) pltMat[numOflevel+1,] <- naCount
  
  # transform to relative frequency
  pltMat <- pltMat %*% diag(1/colSums(pltMat))
  return(pltMat)
}

### Plot global improvment,LBP intensity, Roland-morris Score
# cluster: obs cluster membership
# valData: data set only contains first 2-10 vars
plotValidateVar <- function(valData, cluster){
  
  df <- plotClust(cluster, valData, c(8,5,2), "IMPROVMENT", 1) 
  p1 <- ggplot(df, aes(Variable, IMPROVMENT, group=Class, color= Class)) + geom_point() + geom_line() + theme(legend.direction = 'horizontal', legend.position = 'top')
  df <- plotClust(cluster, valData, c(9,6,3), "LBPINTENSITY", 1)
  
  p2 <-ggplot(df, aes(Variable, LBPINTENSITY, group=Class, color= Class)) + geom_point() + geom_line() + theme(legend.direction = 'horizontal', legend.position = 'top')  
  
  df <- plotClust(cluster, valData, c(10,7,4), "RolandMorrisScore", 1)
  p3 <-ggplot(df, aes(Variable, RolandMorrisScore, group=Class, color= Class)) + geom_point() + geom_line() + theme(legend.direction = 'horizontal', legend.position = 'top')
  grid.arrange(p1, p2, p3, ncol=3)
}

### plot continuous data: age, height, bmi etc
# cluster: obs cluster membership
# noImputeData: orginal data set
# varName: the name of the variable
plotContVar <- function(noImputeData, cluster, varName, realName){
  
  gps <- length(unique(cluster))
  Class <- factor(cluster, labels = paste("Cluster", 1:gps))
  
  varIndx <- which(names(noImputeData) == varName)
  noImputeData <- as.data.frame(cbind(noImputeData,Class))
  
  # plot
  p <- ggplot(noImputeData, aes(Class,noImputeData[,varIndx])) + geom_boxplot(aes(fill = Class))
  p + scale_fill_brewer(palette="GnBu") + ylab(realName)
}
############################ FUNCTIONS END ###############################

############################ DATA READING START ############################
# read the data set and variable description excel files
rawData <- read.xls("data_challenge.xlsx")
varDpt <- read.xls("variables_data_challenge.xlsx", sheet="Variables") 
names(varDpt) <- c("index", "variable", "domain","shortQ","fullQ","content",
                   "type","N","missing","missPcent","fullContent","reference","comment","typeQ")
############################ DATA READING END ##############################


############################ DATA CLEANING START ############################
data <- rawData

#### INPUTE FALSE NA'S START ####
## work related
data$fabq60 <- ifelse(is.na(data$fabq60) & data$barb0 %in% c(4,5,6,7), -1, data$fabq60)
data$fabq70 <- ifelse(is.na(data$fabq70) & data$barb0 %in% c(4,5,6,7), -1, data$fabq70)

data$fabq80 <- ifelse(is.na(data$fabq80) & data$barb0 %in% c(4,5,6,7), -1, data$fabq80)
data$fabq90 <- ifelse(is.na(data$fabq90) & data$barb0 %in% c(4,5,6,7), -1, data$fabq90)
data$fabq100 <- ifelse(is.na(data$fabq100) & data$barb0 %in% c(4,5,6,7), -1, data$fabq100)
data$fabq110 <- ifelse(is.na(data$fabq110) & data$barb0 %in% c(4,5,6,7), -1, data$fabq110)
data$fabq120 <- ifelse(is.na(data$fabq120) & data$barb0 %in% c(4,5,6,7), -1, data$fabq120)
data$fabq130 <- ifelse(is.na(data$fabq130) & data$barb0 %in% c(4,5,6,7), -1, data$fabq130)
data$fabq140 <- ifelse(is.na(data$fabq140) & data$barb0 %in% c(4,5,6,7), -1, data$fabq140)

## narrow down the levels of variables barbo(working status) to 3 levels
## barbo = c(1,2,3)   => 1 working
## barbo = c(4,5,6,7) => 2 no-working
## barbo = c(8)       => 3 others
data$barb0 <- ifelse(data$barb0 %in% c(1,2,3), 1, data$barb0)
data$barb0 <- ifelse(data$barb0 %in% c(4,5,6,7), 2, data$barb0)
data$barb0 <- ifelse(data$barb0 == 8, 3, data$barb0)

# Pain related, domin_hp = 1 mean no dominate back pain
data$facetextrot <- ifelse(is.na(data$facetextrot) & data$domin_bp == 1, -1, data$facetextrot)
data$facetsit <- ifelse(is.na(data$facetsit) & data$domin_bp == 1, -1, data$facetsit)
data$facetwalk <- ifelse(is.na(data$facetwalk) & data$domin_bp == 1, -1, data$facetwalk)
data$paraspin_debut <- ifelse(is.na(data$paraspin_debut) & data$domin_bp == 1, -1, data$paraspin_debut)

## muscle palpation related. musclepalp = 1 means having replation pain during palpation
## We supsect there are relationship between "triggerpoint", "musclegroup_palp" and "musclepalp".
## so we check the correlation among them. It turn out the musclepalp and musclegroup_palp are highly
## corrlated. So we decide to use musclegroup_palp to update the NA's in musclepalp
table(data$triggerpoint, data$musclegroup_palp)
table(data$musclepalp, data$triggerpoint)
table(data$musclepalp, data$musclegroup_palp)
data$musclegroup_palp <- ifelse(is.na(data$musclegroup_palp) & data$triggerpoint == 0 & data$musclepalp == 0, -1, data$musclegroup_palp)
data$musclepalp <- ifelse(is.na(data$musclepalp) & data$musclegroup_palp %in% c(1,2,3), 1, data$musclepalp)

# check the #NA's after updating the False NA's
length(data$fabq60[is.na(data$fabq60)])
length(data$fabq70[is.na(data$fabq70)])
length(data$fabq80[is.na(data$fabq80)])
length(data$fabq90[is.na(data$fabq90)])
length(data$fabq100[is.na(data$fabq100)])
length(data$fabq110[is.na(data$fabq110)])
length(data$fabq120[is.na(data$fabq120)])
length(data$fabq130[is.na(data$fabq130)])
length(data$fabq140[is.na(data$fabq140)])
length(data$facetextrot[is.na(data$facetextrot)])
length(data$facetsit[is.na(data$facetsit)])
length(data$facetwalk[is.na(data$facetwalk)])
length(data$musclegroup_palp[is.na(data$musclegroup_palp)])
#### INPUTE FALSE NA'S END ####

#### TRANSFORM CATEGROICAL VARS TO FACTOR START ####
# treat all the ordinal variables as continuous
catIndx <- which(varDpt$type %in% c("Dichotomous","Multistate nominal", "Trichotomous"))
data[,catIndx] <- lapply(data[,catIndx], as.factor)
#### TRANSFORM CATEGROICAL VARS TO FACTOR END ####

#### IMPUTE NA'S BY MICE PKG START ####
# drop the first 1 to 10 variables
valData <- data[,1:10]
data <- data[,-c(1:10)]

# count the frequency for the new data set(without first 10 variables)
vars <- dim(data)[2]
obs <- dim(data)[1]
naCountCol <- sapply(1:vars, function(index,o) length(which(is.na(o[,index]))), data)
naCountRow <- sapply(1:obs, function(index,o) length(which(is.na(o[index,]))), data)

# drop the variables with over 30% missing values, which is only "triggerpoint"
cutoff <- 0.30
inxVar <- which(naCountCol >= obs*cutoff)
trpoint <- data[,inxVar]
data <- data[, -inxVar]

# drop the observation with over 30% missing value here(13 obs)
inxRow <- which(naCountRow >= vars*cutoff)
dropObs <- data[inxRow, ]
data <- data[-inxRow,]

valDropObs <- valData[inxRow, ]
valData <- valData[-inxRow,] #also drop them from the validate data set

# we reassign the row index since we drop some obs and col
row.names(data) <- 1:915
noImputeData <- data

# use to plot ord as cate vars
ordPlotData <- data
ordIndx <- which(varDpt$type == "Ordinal") - 10
ordPlotData[, ordIndx] <- lapply(data[,ordIndx], as.factor)

# impute NA's via mice package
#mice.out <- mice(data, m=5, maxit=10, seed=12345)
#save(mice.out, file = "miceFullOrdinal.Rdata")
#save(mice.out, file = "micePartialOrdinal.Rdata")

# consolidate the five imputed data set into one
#load("micePartialOrdinal.Rdata")
#data<-consolDataSet(mice.out, data)

#save(data, file = "consolFullOrdinal.Rdata")
#save(data, file = "consolPartialOrdinal.Rdata")
#### IMPUTE NA'S BY MICE PKG END ####
############################ DATA CLEANING END ###############################

########################## CLUSTERING ANALYSIS STARTS #######################
#### Applied MCA starts ####
library(MASS)
library(cluster)
library(FactoMineR)

load("consolFullOrdinal.Rdata")
#load("consolPartialOrdinal.Rdata")

# Separate the continuous and categorical variables
catIndx  <- unlist(lapply(1:111, function(index,o) class(o[,index])=="factor", data))
cateData <- data[,catIndx]  # pure categorical data set
contData <- data[,!catIndx] # pure continuous data set

# dimension of the pure categorical data set
len <- length(cateData)
mca.out <- MCA(cateData, graph = FALSE, ncp = len - 1)

# decide the number of p
m.eig <- mca.out$eig$eigenvalue
m.eig <- ifelse(m.eig > 1/len, ((len/(len-1))*(m.eig-1/len))^2, 0)

plot(1:(len-1), m.eig[1:(len-1)], bty="n", pch=20)
lines(1:(len-1), m.eig[1:(len-1)])

# attach MCA's coordinates to continuous data
pickComp <- which(cumsum(m.eig)/sum(m.eig) >= 0.95)
reducedData <- cbind(mca.out$ind$coord[,1:pickComp[1]], contData)
#### Applied MCA ends ####

####### RUN ALGORITHMS START ######
#### Kmeans Start ####
set.seed(12345)
k.out <- kmeans(scale(reducedData),3, iter.max = 50, nstart = 10)
k.cluster <- k.out$cluster
#### Kmeans End ####
########################## CLUSTERING ANALYSIS ENDS #######################

############### SECECT THE IMPORTMANT VARIABLES STARTS####################
# complete contribution table
ncp <- dim(reducedData)[2]
pca.out <- PCA(scale(reducedData), ncp=ncp)

p.eig <- pca.out$eig$eigenvalue
pctEig <- p.eig/sum(p.eig)
cumsum(p.eig)/sum(p.eig)

# pick the important principal component that help separate the observations(Pick dim 1-7)
plot(1:45, p.eig, pch = 20) # observe the elbow point
lines(1:45, p.eig)
# pair plot the cluster according to the PCs, only PC 1 and 2 help separate the clusters
pairs(pca.out$ind$coord[,1:7], col=k.cluster+2)

# calcuate the contribution of each variable to PC1 and PC2
colCutoff <- 2
varContrib <-  pca.out$var$contrib[,1:colCutoff]/100
#completeTable <- cbind(varContrib, varContrib %*% pctEig[1:colCutoff])
completeTable <- cbind(varContrib, varContrib %*% sort(p.eig[1:colCutoff], decreasing = T))
colnames(completeTable) <- c(paste("PC", 1:colCutoff), "WeightedContribution")
rownames(completeTable)[1:7] <- paste("MCA", 1:7) 
completeTable <- completeTable[order(completeTable[,colCutoff+1], decreasing = T), ]

#p.propEig <- p.eig[1:2]/sum(p.eig)
#completeTable <-rbind(completeTable, c(p.propEig, sum(p.propEig)))
write.csv(x = completeTable, file = "complete contribition table.csv")

### categorical contribution table
col <- c(1,5) 
dimContrib <- mca.out$var$contrib[,col]/100
#wtContrib <- data.frame(dimContrib %*% (m.eig/sum(m.eig))[col])
wtContrib <- data.frame(dimContrib %*% m.eig[col])

# remove the level index for each variables names
name <- rownames(wtContrib)
pos <- regexpr(pattern ="_[^_]*$", name)
name <-substr(name,1,pos-1)

# create weighted contribution table, cut off up to start50
rawTable <- cbind(name, data.frame(dimContrib, wtContrib))
cateTable <- aggregate.data.frame(rawTable[,2:4], by=list(rawTable$name), FUN = sum)
names(cateTable) <- c("Name",paste("MCA", col), "WeightedContribution")
cateTable <- cateTable[order(cateTable[,length(col)+2],decreasing = T), ]

#m.propEig <- m.eig[col]/sum(m.eig)
#eigRow <- data.frame("Eigen Prop", m.propEig[1], m.propEig[2], sum(m.propEig))
#names(eigRow) <- names(cateTable)
#cateTable <- rbind(cateTable, eigRow)
cateTable[,2:4] <- round(cateTable[,2:4], 5)
write.csv(x = cateTable, file = "categorical contribition.csv")

# pick the important variables
impCateVar <- as.character(cateTable[which(cateTable[,4] > 0.0001),1])
impContVar <- rownames(completeTable[which(completeTable[,3] >= 0.1),])
impContVar <- impContVar[-grep("M.A", impContVar)]
impVarList <- c(impContVar,impCateVar)
write.csv(x = impVarList, file = "important variables.csv")

### run throught the whole process again by the reduced data set(consist of only important variabls)
impCateDate <- cateData[,which(names(cateData) %in% impCateVar)]
impContData <- contData[,which(names(contData) %in% impContVar)]

len <- length(impCateDate)
imp.mca.out <- MCA(impCateDate, graph = FALSE, ncp = len - 1)

# decide the number of p
imp.p.eig <- imp.mca.out$eig$eigenvalue
imp.p.eig <- ifelse(imp.p.eig > 1/len, ((len/(len-1))*(m.eig-1/len))^2, 0)

plot(1:(len-1), m.eig[1:(len-1)], bty="n", pch=20)
lines(1:(len-1), m.eig[1:(len-1)])

# attach MCA's coordinates to continuous data
pickComp <- which(cumsum(imp.p.eig)/sum(imp.p.eig) >= 0.95)
impRedData <- cbind(imp.mca.out$ind$coord[,1:pickComp[1]], impContData)

#impRedData <- impContData  
set.seed(12345)
imp.k.out <- kmeans(scale(impRedData),3, iter.max = 50, nstart = 10)
imp.k.cluster <- imp.k.out$cluster

library(MixGHD)
ARI(imp.k.cluster, k.cluster) ## the clustermember are highly matched.

require(factoextra)
fviz_cluster(imp.k.out, scale(impRedData), geom = "point", xlab = F, ylab=F, main="Projection due to Important Variables")
############### SECECT THE IMPORTMANT VARIABLES ENDS  ####################

############## PLOTTING THE VARIABLES STARTS #############################
#### Continuous variables ####
plotContVar(noImputeData, k.cluster, "age", "Age")
plotContVar(noImputeData, k.cluster, "bhoej0", "Height")
plotContVar(noImputeData, k.cluster, "bmi", "BMI")
plotContVar(noImputeData, k.cluster, "vasl0", "LBP Intensity at Consulation")   # Lbp intensity during consulation
plotContVar(noImputeData, k.cluster, "okon0", "Ability to Decrease Pain")   # able to decrese pain
plotContVar(noImputeData, k.cluster, "obeh0", "Is Treatment Essential")   # Treatment not essential
plotContVar(noImputeData, k.cluster, "htil0", "Self-rate General Health")   # Self-rated general health
plotContVar(noImputeData, k.cluster, "rmprop", "Roland-Morris summary score")   # Roland-Morris summary score
#### Continuous variables ####

#### Ordinal variables ####
# pain dlva0
pltMat <- plotCateVar(ordPlotData, "dlva0", k.cluster)
rownames(pltMat) <- c("1-2 weeks", "2-4 weeks", "1-3 months", ">3 months", "Missing")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("duration","class","frequency")
ggplot(pltDf, aes(duration,frequency, group=class, fill=duration)) + geom_bar(stat="identity") + facet_wrap(~ class) + theme(legend.direction = 'horizontal', legend.position = 'top')

# Felt lacking in energy and strength: mid3
plotContVar(noImputeData, k.cluster, "mdi3", "Felt lacking in energy and strength")

pltMat <- plotCateVar(ordPlotData, "mdi3", k.cluster)
rownames(pltMat) <- c("no time", "some time", "< half", "> half","mostTime", "all time", "NA")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("LackinginEnergy","class","frequency")
ggplot(pltDf, aes(LackinginEnergy,frequency, group=class, fill=LackinginEnergy)) + geom_bar(stat="identity") + facet_wrap(~ class)

#less self-confidence Mdi4
plotContVar(noImputeData, k.cluster, "mdi4", "less self-confidence")

pltMat <- plotCateVar(ordPlotData, "mdi4", k.cluster)
rownames(pltMat) <- c("at no time", "some time", "slightly-all the time", "Missing")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("lessSelfConfidence","class","frequency")
ggplot(pltDf, aes(lessSelfConfidence,frequency, group=class, fill=lessSelfConfidence)) + geom_bar(stat="identity") + facet_wrap(~ class)

# fabq120
plotContVar(noImputeData, k.cluster, "fabq120", "Cann't Work with Current pain")

pltMat <- plotCateVar(ordPlotData, "fabq120", k.cluster)
rownames(pltMat) <- c("No Working",paste("Pain", 0:6), "NA")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("UnableToWorkWithPain","class","frequency")
ggplot(pltDf, aes(UnableToWorkWithPain,frequency, group=class, fill=UnableToWorkWithPain)) + geom_bar(stat="identity") + facet_wrap(~ class)

# bfbe0, physicial loading at work, !!!!
plotContVar(noImputeData, k.cluster, "bfbe0", "physical loading at work")

pltMat <- plotCateVar(ordPlotData, "fabq120", k.cluster)
rownames(pltMat) <- c("No Working",paste("Pain", 0:6), "NA")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("UnableToWorkWithPain","class","frequency")
ggplot(pltDf, aes(UnableToWorkWithPain,frequency, group=class, fill=UnableToWorkWithPain)) + geom_bar(stat="identity") + facet_wrap(~ class)
#### Ordinal variables ####

#### Categorical variables ####
# sex
pltMat <- plotCateVar(noImputeData, "bsex0", k.cluster)
rownames(pltMat) <- c("female", "male")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("sex","class","frequency")
ggplot(pltDf, aes(sex,frequency, group=class, fill=sex)) + geom_bar(stat="identity") + facet_wrap(~ class)

# start_risk, summary score of start10~90
pltMat <- plotCateVar(noImputeData, "start_risk", k.cluster)
rownames(pltMat) <- c("low risk", "medium risk", "hight risk", "Missing")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("HightRiskGroup","class","frequency")
ggplot(pltDf, aes(HightRiskGroup, frequency, group=class, fill=HightRiskGroup)) + geom_bar(stat="identity") + facet_wrap(~ class) + theme(legend.direction = 'horizontal', legend.position = 'top')

# romflex: Pain on flexion (AROM)
pltMat <- plotCateVar(noImputeData, "romflex", k.cluster)
rownames(pltMat) <- c("low pain", "back pain", "LPainOrNBackkPain", "Missing")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("PainOnFlexion","class","frequency")
ggplot(pltDf, aes(PainOnFlexion, frequency, group=class, fill=PainOnFlexion)) + geom_bar(stat="identity") + facet_wrap(~ class)

# working status
# we use the raw data to plot working status since we rescale it before
dropIndx <- attributes(dropObs)$row.names   # drop obs's index
workStatus <- rawData[-dropIndx,]           # drop obs's from raw data
workStatus <- workStatus[,-c(1:10)]         # drop the first 10 variables
row.names(workStatus) <- 1:915              # reindex the row id
workStatus[,4] <- as.factor(workStatus[,4]) # make col of working status as factor type

pltMat <- plotCateVar(workStatus, "barb0", k.cluster)
rownames(pltMat) <- c("self emp", "fl-time", "pt-time","stud","unEmp","early rtir", "hlth rtir","other","missing")
colnames(pltMat) <- paste("cluster", 1:length(unique(k.cluster)))

pltDf <- melt(pltMat)
names(pltDf) <- c("status","class","frequency")
ggplot(pltDf, aes(status,frequency, group=class, fill=status)) + geom_bar(stat="identity") + facet_wrap(~ class) + theme(legend.direction = 'horizontal', legend.position = 'top')
#### Categorical variables ####

#### Validated variables ####
plotValidateVar(valData,  k.cluster)
#### Validated variables ####

#### plot important variables based on the cluster membership ####
varIndx <- which(names(noImputeData) %in% impContVar)
graphName <- "ClusterDifference"
tempData <- noImputeData
tempData[,varIndx] <- scale(tempData[,varIndx])
df <- plotClust(k.cluster,tempData,varIndx, graphName, 2)
p1 <- ggplot(df, aes(Variable, ClusterDifference, group=Class, color= Class))
p1 + geom_point() + geom_line() + ylab("Cluster Difference")
#### plot important variables based on the cluster membership ####

#### projection ####
# choose.vars = c("Dim 1", "Dim 5", impContVar)
fviz_cluster(k.out, scale(reducedData), geom = "point")

fviz_pca_biplot(
  prcomp(scale(reducedData)), axes=c(1,2),
  habillage = k.cluster,
  addEllipses = F,
  col.var = "grey",
  #alpha.var = "cos2",
  label = "var",
  title=""
) + scale_color_brewer(palette = "Dark2") + theme_minimal() + theme(legend.direction = 'horizontal', legend.position = 'top')
#### projection ####
############## PLOTTING THE VARIABLES END #############################

############## OUTPUT THE CLUSTER MEMBERSHIP STARTS ############################
output <- data.frame(cbind(rawData[,1], rep("Excluded due to too many missing values",length=928)), stringsAsFactors=F)
output[-dropIndx,2] <- paste("Cluster", k.cluster)
colnames(output) <- c("ID", "Cluster Membership")
write.csv(output, file = "cluster membership.csv")
############## OUTPUT THE CLUSTER MEMBERSHIP ENDS ########################