
debug <- FALSE
if ( debug ) {
    library("Rcpp")
    source("~/programs/segmenTier/pkg/R/plot.R")
    source("~/programs/segmenTier/pkg/R/cluster.R")
    source("~/programs/segmenTier/pkg/R/segment.R")
    sourceCpp("~/programs/segmenTier/pkg/src/segment.cpp")
    sourceCpp("~/programs/segmenTier/pkg/src/cluster.cpp")
    load("~/programs/segmenTier/pkg/data/primseg436.rda")
} else {

    ## load library
    library("segmenTier")
    ## load time-series data
    ## contains tsd from primseg436 for
    ## a 7.6 kb genomic region
    data(primseg436)
}

plot2file <- !interactive() # plot to pdf if called from command-line
fig.type <- "eps"

## EXAMPLE DATASET FROM BUDDING YEAST
## NOTE Here we vary the parameters with strongest effect
## on segmentation of our budding yeast data set.
## Parameters with little effect on similarity-based scoring functions
## icor and ccor are discussed in demo(segment_test)


### TIME-SERIES PROCESSING PARAMETERS
## NOTE that the current pipe-line for batch processing
## allows only one configuration of time-series processing
## while all downstream steps (clustering, segmentation)
## can be varied over parameter ranges.
## NOTE that currently segmenTier is not tested for clustering
## of time-series directly (without Discrete Fourier Transform)
## but can be tested by setting use.fft to FALSE
trafo <- "raw"     # transformation of the raw data 
use.fft <- TRUE    # cluster discrete Fourier transform of data?
use.snr <- TRUE    # use DFT scaling (SNR is described as
                   # relative amplitude scaling in Machne&Murray, PLoS ONE 2012)
dft.range <- 1:7   # range of DFT to use for clustering
dc.trafo <- "ash"  # transformation of the first (DC) component of the DFT
                   # NOTE: add component 1 (DC) to DFT range to use
low.thresh <- -Inf # minimal total signal (DC component of DFT if use.fft)

### CLUSTERING PARAMETERS
K <- c(12)         # cluster number K; multiple allowed; specifically, note
                   # that k-means has a random effect at initialization
                   # and replicates of the same K can  yield different
                   # results for otherwise 
nui.thresh <- 0.6  # threshold of position-cluster correlation below which
                   # the position will be assigned to the nuisance cluster
## k-means initialization
iter.max <- 100000 # max. iterations in kmeans
nstart <- 100      # number of initial configurations tested in kmeans

### SEGMENTATION PARAMETERS
## segmenTier parameters are handled via the settings function,
## where all parameters can be passed as vectors.
vary <- setVarySettings(
    E=c(1,3), # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=c(150), # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=c(1,3)#-/+ correlation of nuisance cluster with others and itself
)

## PRE-PROCESS TIME SERIES FOR CLUSTERING
## take DFT and scale amplitudes, and
## select components of DFT
tset <- processTimeseries(ts=tsd, na2zero=TRUE,
                          trafo=trafo, dc.trafo=dc.trafo,
                          use.fft=use.fft, dft.range=dft.range,
                          use.snr=use.snr, low.thresh=low.thresh)

## CLUSTER PRE-PROCESSED TIME SERIES
set.seed(15) # stable kmeans clustering
cset <- clusterTimeseries(tset, K=K, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)


### CALCULATE SEGMENTS FOR ALL CLUSTERINGS and
### FOR CHOSEN SEGMENTATION PARAMETERS
## NOTE, that the function optionally also stores the
## the scoring and backtracing matrices (see demo/segment_test.R)
sset <- segmentCluster.batch(cset, varySettings=vary,
                             id="mysegments",
                             type.name=c("E","M","nui"), # segment type names
                             verb=1, save.matrix=FALSE) 

## INSPECT RESULTS

## ... the print method shows the segment borders and cluster associations
print(sset)

## the segment coordinates are found in:
head(sset$segments)

## plot segmentation
if ( plot2file )
  plotdev("segment_data",res=300,width=10,height=5,type=fig.type)

# plot.matrix=TRUE will additionally plot the internal scoring matrices
# if segmentation was calculated with save.matrix=TRUE
plotSegmentation(tset, cset, sset, cex=.5, lwd=2) 

if ( plot2file )
  dev.off()



## BEST AND WORST PARAMETER SETS, as resulting from parameter scan
## in publication
## UNDER-FRAGMENTATION - red cluster in fig 2
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=200,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=1   #-/+ correlation of nuisance cluster with others and itself
)
bad1 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))
## OVER-FRAGMENTATION - magenta
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=3   #-/+ correlation of nuisance cluster with others and itself
)
bad2 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## as above but with ccor
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="ccor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=3   #-/+ correlation of nuisance cluster with others and itself
)
bad2.ccor <- segmentCluster.batch(cset, varySettings=vary,type.name=c("S"))

## BEST FRAGMENTATION - cyan

## best 1 - high E/nui - long M
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=200,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=3   #-/+ correlation of nuisance cluster with others and itself
)
best1 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))
## best 2 - intermediate
vary <- setVarySettings(
    E=2,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=3   #-/+ correlation of nuisance cluster with others and itself
)
best2 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))
## best 3 - E/nui=1 - short M
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=1   #-/+ correlation of nuisance cluster with others and itself
)
best3 <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))
## best 3-ccor - E/nui=1 - short M
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="ccor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=1   #-/+ correlation of nuisance cluster with others and itself
)
best3.ccor <- segmentCluster.batch(cset, varySettings=vary,type.name=c("S"))

## use layout to combine plots
#plotdev("segment_data_examples",res=300,width=10,height=5,type="jpeg")
if ( plot2file ) # Figure 3 of the preprint manuscript
  plotdev("segment_data_examples",res=300,width=10,height=5,type=fig.type)
layout(matrix(1:10,ncol=1),heights=c(.25,.5,.5,.075,.075,.075,.075,.075,.075,.075))
par(mai=c(0.1,2,0.05,0.01),xaxs="i",yaxs="r")
par(cex=1) 
plot(tset,ylabh=TRUE)
par(cex=.6) 
plot(cset,axes=2,cex=.7,ylabh=FALSE); mtext("clustering",2,7.2,las=2,cex=1.2)
par(cex=1.2) # increase axis labels
par(mai=c(0.01,2,0.01,0.01))
plot(bad1,"segments",lwd=3)
plot(best3,"segments",lwd=3)
plot(best3.ccor,"segments",lwd=3)
plot(best2,"segments",lwd=3)
plot(best1,"segments",lwd=3)
plot(bad2,"segments",lwd=3)
plot(bad2.ccor,"segments",lwd=3)
if ( plot2file )
  dev.off()

## SYSTEMATIC VARIATION OF SPECIFIC PARAMETERS

## NOTE: the following produces Figure S4a of the preprint manuscript


## vary M; E=nui=2
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=seq(50,250,25), # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=1   #-/+ correlation of nuisance cluster with others and itself
)
varM <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## vary E; nui=1, M=150
vary <- setVarySettings(
    E=1:9,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=2  #-/+ correlation of nuisance cluster with others and itself
)
varE <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## vary nui; E=1, M=150
vary <- setVarySettings(
    E=1,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=150,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=1:9   #-/+ correlation of nuisance cluster with others and itself
)
varN <- segmentCluster.batch(cset, varySettings=vary,type.name=c("E","M","nui"))

## NOTE: use layout to combine plots
if ( plot2file ) # Figure S4a of the preprint manuscript
  plotdev("segment_data_scans",res=300,width=10,height=7.5,type=fig.type)
layout(matrix(1:4,ncol=1),heights=c(.5,.5,.5,.5))
par(mai=c(0.1,2,0.05,0.01),xaxs="i",yaxs="r")
#par(cex=1) 
#plot(tset,ylabh=TRUE)
par(cex=.6) 
plot(cset,axes=2,cex=.7)
par(cex=1.2) # increase axis labels
par(mai=c(0.0,2,0.0,0.01))
plot(varM,"segments",lwd=3)
plot(varE,"segments",lwd=3)
plot(varN,"segments",lwd=3)
if ( plot2file )
  dev.off()

### MULTIPLE CLUSTERINGS
## Here we generate multiple clusterings, and segmentations
## (here with constant parameters) will be calculated for all of them.
## NOTE that nui.thresh acts as a noise filter for clustering, 
## based on a minimal position-cluster similarity in matrix cset$Pci
## NOTE that random effects of k-means clustering could potentially
## be utilized to clean data.
## cluster
nui.thresh <- nui.thresh
vK <- rep(12,6) # c(16,16,16) #,20,20,20)
set.seed(10) # stable kmeans clustering
kset <- clusterTimeseries(tset, K=vK, iter.max=iter.max, nstart=nstart,
                          nui.thresh=nui.thresh)

## calculate segments
vary <- setVarySettings(
    E=3,    # scale exponent of similarity matrices csim
    S="icor", # SCORING FUNCTIONS
    M=75,   # scoring function minimal length penalty
    Mn=100,   # M for nuisance clusters
    nui=3   #-/+ correlation of nuisance cluster with others and itself
)
#vary$nui <- vary$E <- 3
vark <- segmentCluster.batch(kset, varySettings=vary)

## plot segmentations
if ( plot2file ) # Figure S4b of the preprint manuscript
    plotdev("segment_data_clusterings",res=300,width=10,height=7.5,type=fig.type)
plotSegmentation(NULL, kset, vark, cex=.5, lwd=2, mai=c(0.1,2,0.05,0.01)) 

if ( plot2file )
    dev.off()


## TODO: PROOF OF CONCEPT: CLUSTER-FREE SEGMENTATION

## CLUSTER-FREE - NOTE: this takes a lot of memory and time!
## each position x_i is treated as its own cluster
##N <- nrow(tset$dat)
##D <- sum(!tset$rm.vals)
##seq <- rep(0, N)
##seq[!tset$rm.vals] <- 1:D
##P <- matrix(NA,nrow=N,ncol=D)
#### position-position cross-correlation
##P[!tset$rm.vals,] <- clusterCor_c(tset$dat[!tset$rm.vals,],
##                                  tset$dat[!tset$rm.vals,])
##cset <- list()
##cset$clusters <- matrix(seq,ncol=1)
##cset$Pci <- list(P)
##colnames(cset$clusters) <- names(cset$Pci) <- "all"
##
##vary$S <- "icor" # only 'icor' is possible for this approach!
##vary$M <- 150
##vary$nui <- vary$E <- 3
##sset <- segmentCluster.batch(cset, varySettings=vary)
