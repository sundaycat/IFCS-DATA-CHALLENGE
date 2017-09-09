
###########################################################################
####### DETERMINE OPTIMAL NUMBER OF CLUSTERS G - BEGIN ####################
###########################################################################

# run clustering algorithms to figure out optimal G:  gpcm, msn, GHD, 
library(mixture)
library(EMMIXskew)
library(MixGHD)

gmm.out=gpcm(as.matrix(reducedData),G=3:12,mnames="VVV")
gmm.out$BIC  #pick G that minimizes BIC

# library(teigen)
# set.seed(12345)
# reducedData must be scaled
# tdist=teigen(x=scale(reducedData),Gs=3:12,models="UUUU", maxit = c(30,30)) 
# NOTE:  teigen did not converge for G=3:12

set.seed(12345)
BIC <- NULL
msn.out <- list() # set msn.out to be a list
for(g in 3:8){
  msn.out[[g]] <- EmSkew(dat=reducedData,g=g,distr="msn",debug=FALSE)  # stores the output from EmSkew in a list designated by the number of groups
  BIC[g] <- msn.out[[g]]$bic
}
# NOTE:  this produce error when G=3:12 but works when G=3:5
# model suggests G=3

set.seed(12345)
GHD.out=MGHD(reducedData,G=3:6,modelSel="BIC") # model suggests G=3

set.seed(12345)
MCGHD.out=MCGHD(reducedData,G=3:6,modelSel="BIC") # model suggests G=3

set.seed(12345)
MSGHD.out=MSGHD(reducedData,G=3:6,modelSel="BIC") # model suggests G=3

##############################################################################
####### DETERMINE OPTIMAL NUMBER OF CLUSTERS G - END #########################
##############################################################################

##############################################################################
#### COMPARING CLUSTER SOLUTIONS USING KNOWN NUMBER OF CLUSTERS - BEGIN ######
##############################################################################

# attach MCA's coordinates to continuous reducedData
reducedData <- cbind(mca.out$rs[,1:11], contreducedData)

# run several clustering algorithms to compare
G <- 3
library(cluster)
library(fclust)
library(FPDclustering)
library(EMMIXskew)
library(mixture)
library(MixGHD)
library(mmtfa)
library(FisherEM)

# run clsutering algorithms:  k-means, pam, fkm, pdclust, mvn, mst, GHD, MTFA, DLM
set.seed(12345)
kmeans.out <- kmeans(reducedData, centers=G, iter.max = 10, nstart=10)
pam.out <- pam(reducedData, k=G)
fkm.out <- FKM(reducedData, k=G)
pdc.out <- PDclust(reducedData, k=G)

mvn.out <- EmSkew(reducedData, G, distr = "mvn", debug = FALSE)
mst.out <- EmSkew(reducedData, G, distr = "mst", debug = FALSE)

GHD.out <- MGHD(reducedData, G=G, modelSel = "BIC")

MTFA.out <- mmtfa(reducedData, models="UUUU", Gs=G, init="hard")

DLM.out <- fem(reducedData, G)

# compare clutering solutions using ARI
ARI(kmeans.out$cluster, pam.out$clustering)        # 0.806020
ARI(kmeans.out$cluster, fkm.out$clus[,1])          # 0.8735205
ARI(kmeans.out$cluster, pdc.out$label)             # 0.5292999
ARI(kmeans.out$cluster, mvn.out$clust)             # 0.03591541

ARI(pam.out$clustering, fkm.out$clus[,1])          # 0.7077434
ARI(pam.out$clustering, pdc.out$label)             # 0.5855523

ARI(fkm.out$clus[,1], pdc.out$label)               # 0.5403499

ARI(mvn.out$clust, mst.out$clust)                  # 0.7129204
ARI(mvn.out$clust, GHD.out@map)                    # 0.5147287

ARI(mst.out$clust, GHD.out@map)                    # 0.683744

# compare clustering solutions with Silhouette
sval.kmeans <- silhouette(x=kmeans.out$cluster, dist=dist(reducedData))
plot(sval.kmeans, col="blue") # avg width = 0.29
sval.fkm <- silhouette(x=fkm.out$clus[,1], dist=dist(reducedData))
plot(sval.fkm, col="blue") # avg width = 0.28
sval.pam <- silhouette(x=pam.out$clustering, dist=dist(reducedData))
plot(sval.pam, col="blue")
sval.pdc <- silhouette(x=pdc.out$label, dist=dist(reducedData))
plot(sval.pdc) # avg width = 0.16

##############################################################################
#### COMPARING CLUSTER SOLUTIONS USING KNOWN NUMBER OF CLUSTERS - END ########
##############################################################################