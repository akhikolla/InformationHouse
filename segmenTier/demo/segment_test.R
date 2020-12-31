
#library("segmenTier")

debug <- FALSE#TRUE#
if ( debug ) {
    library("Rcpp")
    source("~/programs/segmenTier/pkg/R/plot.R")
    source("~/programs/segmenTier/pkg/R/cluster.R")
    source("~/programs/segmenTier/pkg/R/segment.R")
    sourceCpp("~/programs/segmenTier/pkg/src/segment.cpp")
} else {
    library("segmenTier")
}

## A sequence of clusters: this must be of type "numeric", and
## most importantly, a nuisance cluster `0' indicates
## noisy or absent data. The nuisance cluster is handled
## with separate parameters, which mostly helps to define
## "tighter" ends of segments with similarity-based scoring functions!

seq <- c(4,4,1,2,4,4,4,2,4,4,3,4,4,2,4,1,4,4,0,0,0,0,0,0,1,2,2,2,4,1,1,1,1,
         0,0,1,1,1,1,1,3,3,3,0,3,3,3,3,0,1,3,4,3,2,4,4,1,3)


## SCORING FUNCTION "ccls": score only by cluster membership
## this is all we need for a segmentation with the simplest
## scoring function "ccls" which is defined by three parameters
sset <- segmentClusters(seq = seq,
                        S = "ccls", M = 3, Mn = 3, a = -2, 
                        save.matrix = TRUE, rm.nui= FALSE)
## the returned structure has class "segments"
class(sset)

## ... for which a plot method is defined that can plot the segments, and,
## if option save.matrix was set to TRUE,
## the internal scoring matrices `S1(i,c)` and `S(i,c)`
par(mfcol=c(3,1),mai=c(0,2,0,0))
plot(sset, plot=c("S1","S", "segments"), lwd=3)

## ... the print method shows the segment borders and cluster associations
print(sset)

## the segment coordinates are found in:
head(sset$segments)


## CLUSTER SIMILARITIES

## For illustration, we manually create a "clustering" object
## as returned by function clusterTimeseries later. It mainly
## comprises a matrix of one or more clusterings and for
## each clustering, the cluster similarity matrices required
## by scoring functions ccor and icor.

## get list of clusters
C <- sort(unique(seq))
C <- C[as.character(C)!="0"]

## SCORING FUNCTION "ccor": cluster-cluster similarity (correlation)
## list of non-nuisance clusters
Ccc <- matrix(0,ncol=length(C),nrow=length(C))
colnames(Ccc) <- rownames(Ccc) <- as.character(C)
diag(Ccc)<- 1 # set diagonal to 1
Ccc[1,2] <- Ccc[2,1] <- -.6 # just enough to avoid merging of 2 with 1
Ccc[2,4] <- Ccc[4,2] <- -.4 # just enough to avoid merging of 2 with 4
Ccc[1,3] <- Ccc[3,1] <- -.5

## SCORING FUNCTION "icor": position-cluster similarity (correlation)
Pci <- matrix(0,ncol=length(C),nrow=length(seq))
set.seed(42) # to keep results constant
for ( i in 1:length(seq) ) {
    for ( j in 1:ncol(Pci) )
        Pci[i,j] <- sample(seq(-1,.1,.01))[3] # set bad correlation
    if ( seq[i]!=0 ) { # set good correlation to its own cluster
        Pci[i,seq[i]] <- sample(seq(.3,1,.01))[1]
    }
}

## construct "clustering" set manually:
cset <- list()
class(cset) <- "clustering"
cset$clusters <- matrix(seq,ncol=1) # a matrix of one or more clusterings
colnames(cset$clusters) <- paste("K",length(C),sep="") # an ID

cset$Ccc <- cset$Pci <- list() # similarity matrices for scoring functions
cset$Ccc[[1]] <- Ccc # ccor: cluster-cluster similarity
cset$Pci[[1]] <- Pci # icor: position-cluster similarity
names(cset$Ccc) <- names(cset$Pci) <- colnames(cset$clusters)

## CLUSTER SORTING & COLORING
## add sorting and coloring to "clustering" object
## clusters are sorted sequentially via similarity matrix Ccc,
## see ?colorClusters
## TODO: align sorting between cset and sset! if cset is unsorted
## 
cset <- colorClusters(cset)
class(cset)
plot(cset) # plot method for class "clustering"

## SIMILARITY BASED SCORING FUNCTIONS
## ccor requires matrix Ccc: cluster-cluster similarity, here correlation
## icor required matrix Pci: position-cluster similarity, here correlation
sset <- segmentClusters(seq = cset,
                        S = "ccor", M = 3, Mn = 3, 
                        save.matrix = TRUE, rm.nui= FALSE)

## Note, that here we keep nuisance segments for illustration.
## Nuisance segments are plotted in gray, while data-based
## segments are colored.

## PLOT FUNCTION FOR CLASS "clustering"
par(mfcol=c(3,1),mai=c(0,2,0,0))
par(xaxs="i") # required to align x-axes with the heatmap plots
cs <- plot(cset) 
plot(sset, plot=c("S", "segments"), lwd=3) # plot segmentation
axis(1)

## PARAMETER SCAN: the batch function

## get settings structure, with all parameters and their defaults
## parameters you want to test can be supplied as numeric or
## character vectors
parameters <- setVarySettings(M=3, Mn=3, a=-2, nui=1, E=1,
                              S=c("ccls","ccor","icor"), #"ccor",#
                              multi=c("max","min"),
                              multib=c("max","min","skip"),
                              nextmax=c(TRUE,FALSE))


sset <- segmentCluster.batch(cset = cset, id="test",
                             varySettings=parameters,
                             save.matrix = TRUE, rm.nui= FALSE)


## PLOT ALL
## for paper as pdf
#pdf("segment_data.pdf", width=10,height=4)
layout(matrix(1:4), heights=c(.15,.25,.25,.25))
par(mai=c(.01,2.5,.05,.1),xaxs="i",yaxs="i")
plot(cset,ylim=c(0.5,4.5))
axis(1,mgp=c(1,.2,0))
for ( S in c("ccls", "ccor", "icor") ) {
    plot(sset, plot="segments", params=c(S=S), lwd=2)
    axis(1,mgp=c(1,.2,0))
}
#dev.off()

## NOTE that the parameters "multi", "multib" and "nextmax" only
## affect the basic scoring function "ccls" which is not recommended
## for real data sets, unless the clustering is well defined, and
## and no continuous cluster similarities can be defined.
## Parameters with effects on real data are discussed in demo(segment_data)


