#' segmenTier : cluster-based segmentation
#' from a sequential clustering
#'@author Rainer Machne \email{raim@tbi.univie.ac.at}, Douglas B. Murray,
#' Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#'@docType package
#'@name segmenTier
#'@section Dependencies: The package strictly depends only on
#' \code{Rcpp}.
#' All other dependencies are usually present in a
#' basic installation (\code{stats}, \code{graphics}, \code{grDevices})).
#' @references
#' Machne, Murray & Stadler (2017) <doi:10.1038/s41598-017-12401-8>,
#' Machne & Murray (2012) <doi:10.1371/journal.pone.0037906>, and
#' Lehmann et al. (2013) <doi:10.1186/1471-2105-14-133>
#'@importFrom Rcpp evalCpp
#'@importFrom stats qt sd var BIC AIC
#'@importFrom graphics image axis par plot matplot points lines legend arrows strheight strwidth text mtext abline polygon
#'@importFrom grDevices png dev.off rainbow gray xy.coords
#'@useDynLib segmenTier
NULL # this just ends the global package documentation



### DYNAMIC PROGRAMMING BASED SEGMENTATION OF A CLUSTERING
### implemented by Rainer Machne, hopefully
### as conceived by Peter F. Stadler

### NOTE: these functions are tested in $GENBRO/src/segment_test.R !

## PROBLEM
## find optimal segments k:i in a sequence of clusters
##      k      i     N
## 1333533335134424413   
## SCORING (dyn.prog.) - dynamically fill the matrix S(i,c):
## for i=1..N   // position
##   for c in 1:4 // clusters 1..5
##     // find the segmentation points k at the maximal sum of
##     // S(k-1,c')    , which is the end of the best prev. local segment c'!=c
##     //   PLUS
##     // score(k,i,c) , which is a measure of cluster c enrichment in k..i
##     S[i,c] <- max_c' max_k S(k-1,c') + score(k,i,c)
##     // and store the used k for back-tracing
##     K[i,c] <- k  // the k which delivered S[i,c]
## BACKTRACING: 
## ... and find segments by back-tracing the cluster c and position
## k which had delivered the maximal score S(c,i) at each i
## find max score in S(i,c) at i=N
## i = N
## while(i != 0)
##   for all k<i and color c'!= c:
##      if S(i,c) == S(k-1,c') + score[k,c,i] )
##      break;
##   output:  interval k to i with color c
##   i <- k - 1 // nextmax: search next non-decreasing S(i,c)
##   c <- c'


### FUNCTIONS

### MESSAGE UTILS

## nicer time-stamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=stdout()) # until piping is implemented
## stored "warnings" (actually detailed messages)
warn <- function(w, warnings,verb=FALSE) {
  if (verb) cat(w)
  c(warnings,w)
}

### HIGH-LEVEL WRAPPERS - TODO

## TODO: high-level wrapper that takes a time-series as input
## and clusters the time-series  calling segmentClusters
## additionally reports data medians/centers for all segments
## and cluster centers; can be used as input for clusterSegments
## to finally cluster all segments genome-wide
segmentData <- function() {}

## TODO: high-level wrapper that segments genome data into primary domains,
## and sub-divides these into coherent segments based on clustering
## of the data and a dynamic programming algo; ....
clusterSegments <- function() {}


### SEGMENTATION BY DYNAMIC PROGRAMMING - MAIN FUNCTIONS
## NOTE that most of the work is done in segment.cpp
## using the Rcpp interface to C++

#' Run the \code{segmenTier} algorithm.
#' 
#' segmenTier's main wrapper interface, calculates segments from a
#' clustering sequence. This will run the segmentation algorithm once
#' for the indicated parameters. The function
#' \code{\link{segmentCluster.batch}} allows for multiple runs over
#' different parameters or input-clusterings.
#'
#' @details This is the main R wrapper function for the `segmenTier'
#' segmentation algorithm. It takes an ordered sequence of cluster
#' labels and returns segments of consistent clusterings, where
#' cluster-cluster or cluster-position similarities are
#' maximal. Its main input (argument \code{seq}) is either a
#' "clustering" object returned by \code{\link{clusterTimeseries}}
#' (scenario I), or an integer vector of cluster labels (scenario
#' II) or. The function then runs the dynamic programming algorithm
#' (\code{\link{calculateScore}}) for a selected scoring function
#' and an according cluster similarity matrix, followed by the
#' back-tracing step (\code{\link{backtrace}}) to find segment
#' borders.
#'
#' The main result, list item "segments" of the returned 
#' object, is a 3-column matrix, where column 1 is the cluster
#' assignment and columns 2 and 3 are start and end indices of the
#' segments. For the batch function \code{\link{segmentCluster.batch}},
#' the "segments" item is a \code{\link[base:data.frame]{data.frame}}
#' contain additional information, see ?segmentCluster.batch.
#' 
#' As shown in the publication, the parameters \code{M},
#' \code{E} and \code{nui} have the strongest impact on resulting
#' segment borders.  Other parameters can be fine-tuned but had
#' little impact on our test data set.
#' 
#' In the default and tested scenario I, when the input is an object
#' of class "clustering" produced by \code{\link{clusterTimeseries}},
#' the cluster-cluster and cluster-position similarity matrices are
#' already provided by this object.
#'
#' In the second scenario II for custom use, argument \code{seq} can
#' be a simple clustering vector, where a nuisance cluster must be
#' indicated by cluster label "0" (zero). The cluster-cluster or
#' cluster-position similarities MUST be provided (argument
#' \code{csim}) for scoring functions "ccor" and "icor",
#' respectively. For the simplest scoring function "ccls", a uniform
#' cluster similarity matrix is constructed from arguments \code{a}
#' and \code{nui}, with cluster self-similarities of 1,
#' "dissimilarities" between different clusters using argument
#' \code{a<0}, and nuisance cluster self-similarity of \code{-a}.
#'
#' The function returns a list (class "segments") comprising of the
#' main result (list item "segments"), and "warnings" from the dynamic
#' programming and backtracing phases, the used similarity matrix
#' \code{csim}, extended by the nuisance cluster; and optionally (see
#' option \code{save.matrix}) the scoring vectors \code{S1(i,c)}, the
#' total score matrix \code{S(i,c)} and the backtracing matrix
#' \code{K(i,c)} for analysis of algorithm performance for novel data
#' sets.  Additional convenience data is reported, such as cluster
#' colors and sortings if argument \code{seq} was of class
#' 'clustering'. These allow for convenient inspection of all data
#' processing steps with the plot methods. A plot method exists that
#' allows to plot segments aligned to "timeseries" and "clustering"
#' plots.
#' @param seq Either an integer vector of cluster labels, or a
#'     structure of class 'clustering' as returned by
#'     \code{\link{clusterTimeseries}}. The only strict requirement
#'     for the first option is that nuisance clusters (which will be
#'     treated specially during the dynamic programming routine) have
#'     to be '0' (zero).
#' @param k if argument \code{seq} is of class 'clustering' the kth
#'     clustering will be used; defaults to 1
#' @param csim The cluster-cluster or position-cluster similarity
#'     matrix for scoring functions "ccor" and "icor" (option
#'     \code{S}), respectively. If \code{seq} is of class 'clustering'
#'     \code{csim} is optional and will override the similarity
#'     matrices in \code{seq}. If argument \code{seq} is a simple
#'     vector of cluster labels and the scoring function is "icor" or
#'     "ccor", an appropriate matrix \code{csim} MUST be
#'     provided. Finally, for scoring function "ccls" the argument
#'     \code{csim} will be ignored and the matrix is instead
#'     automatically constructed from argument \code{a}, and using
#'     argument \code{nui} for the nuisance cluster.
#' @param E exponent to scale similarity matrices
#' @param S the scoring function to be used: "ccor", "icor" or "ccls"
#' @param M segment length penalty. Note, that this is not a strict
#'     cut-off but defined as a penalty that must be "overcome" by
#'     good score.
#' @param Mn segment length penalty for nuisance cluster. Mn<M will
#'     allow shorter distances between "real" segments; only used in
#'     scoring functions "ccor" and "icor"
#' @param a a cluster "dissimilarity" only used for pure cluster-based
#'     scoring w/o cluster similarity measures in scoring function
#'     "ccls".
#' @param nui the similarity score to be used for nuisance clusters
#'     in the cluster similarity matrices
#' @param nextmax go backwards while score is increasing before
#'     opening a new segment, default is TRUE
#' @param multi handling of multiple k with max. score in forward
#'     phase, either "min" (default) or "max"
#' @param multib handling of multiple k with max. score in back-trace
#'     phase, either "min" (default), "max" or "skip"
#' @param rm.nui remove nuisance cluster segments from final results
#' @param save.matrix store the total score matrix \code{S(i,c)} and
#'     the backtracing matrix \code{K(i,c)}; useful in testing stage
#'     or for debugging or illustration of the algorithm;
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @return Returns a list (class "segments") containing the main
#'     result (list item "segments"), and additional information (see
#'     `Details'). A plot method exists that allows to plot clusters
#'     aligned to time-series and segmentation plots.
#' @references Machne, Murray & Stadler (2017)
#'     <doi:10.1038/s41598-017-12401-8>
#' @examples
#' # load example data, an RNA-seq time-series data from a short genomic region
#' # of budding yeast
#' data(primseg436)
#' 
#' # 1) Fourier-transform time series:
#' ## NOTE: reducing official example data set to stay within 
#' ## CRAN example timing restrictions with segmentation below
#' tset <- processTimeseries(ts=tsd[2500:6500,], na2zero=TRUE, use.fft=TRUE,
#'                           dft.range=1:7, dc.trafo="ash", use.snr=TRUE)
#' 
#' # 2) cluster time-series into K=12 clusters:
#' cset <- clusterTimeseries(tset, K=12)
#' 
#' # 3) ... segment it; this takes a few seconds:
#' segments <- segmentClusters(seq=cset, M=100, E=2, nui=3, S="icor")
#' 
#' # 4) inspect results:
#' print(segments)
#' plotSegmentation(tset, cset, segments, cex=.5, lwd=3)
#' 
#' # 5) and get segment border table for further processing:
#' sgtable <- segments$segments
#' 
#' @export
segmentClusters <- function(seq, k=1, csim, E=1,
                            S="ccor", M=175, Mn=20, a=-2, nui=1,
                            nextmax=TRUE, multi="max",multib="max", 
                            rm.nui=TRUE, save.matrix=FALSE, verb=1) {

    ## timing currently only used in verbose mode
    ## TODO: report in results
    ##if ( verb>0 )
    ## start time
    stime <- Sys.time()
    
    ## input: cluster set from clusterTimeseries
    cset <- NULL
    if ( class(seq)=="clustering" ) {
        cset <- seq
        seq <- cset$clusters[,k]
        ## if cset & csim are provided, csim overrides
        ## the cset internal matrix!
        if ( S=="ccor" & missing(csim) ) csim <- cset$Ccc[[k]]
        if ( S=="icor" & missing(csim) ) csim <- cset$Pci[[k]]
    }
    ## input: processed time-series from processTimeseries
    ## TODO: cluster-free approach; generate Pci for all vs. all
           
    ## 1: set up sequence and data
    N <- length(seq)
    seqr <- seq

    ## 1a: map to internal 1:K clustering:
    ## TODO: allow character clusters!?
    map <- sort(unique(seqr)) # clusters
    map <- map[map!=0]        #  nuisance cluster
    names(map) <- map
    map[] <- seq_along(map)
  
    seqr <- map[as.character(seq)]
    seqr[seq==0] <- 0      # original nuisance clusters
    seqr[is.na(seqr)] <- 0 # replace NA by nuisance
        
    ## set-up similarity matrix for ccls
    ## internally 'ccor' is used, and we set up the
    ## cluster-cluster similarity function (matrix) here
    if ( S=="ccls" ) {
        L <- length(unique(seqr))
        csim <- matrix(a, nrow=L, ncol=L) # Delta(C,D!=C) = a
        diag(csim) <- 1 # Delta(C,C) = 1
        nui <- -a # csim will be expanded to contain nuisance below
    }
    
    ## 1b: add nuisance cluster if present:
    ## columns and rows are be added to the similarity matrices, using
    ## nui and -nui as "correlations";
    ## cluster index increased by +1, and the nuisance cluster will be "1"
    ## throughout further processing!

    nui.present <- FALSE
    if ( 0 %in% seqr ) {
        
        nui.present <- TRUE
        ## increase clustering by +1:
        ## nuisance cluster will internally be cluster 1 !
        seqr <- seqr + 1 ## TODO: get rid of this, avoid correction in .cpp
        
        if ( S=="icor" ) {
            ## cor(i,c) - similarity of position i to cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[,as.numeric(names(map))]
            ## add nuisance cluster
            csim <- cbind(rep(-nui,N),csim)
            csim[seqr==1,] <- -nui
            csim[seqr==1,1] <- nui
        }
        if ( S %in% c("ccor") ) {
            ## cor(c,c) - similarity of cluster medians
            ## reduce passed matrix to actually present clusters
            csim <- csim[as.numeric(names(map)),as.numeric(names(map)),
                         drop=FALSE] # allow single cluster? or abort here?
            ## add nuisance cluster
            csim <- rbind(rep(-nui,nrow(csim)+1),
                          cbind(rep(-nui,nrow(csim)), csim))
            csim[1,1] <- nui
        }
        map <- c('0'=1, map  + 1)
    }
    
    ## get clusters
    C <- sort(unique(seqr))

    ## scale similarity matrix!
    sgn <- sign(csim) # store sign
    csim <- csim^E # scale matrix
    ## if exponent is even (checking within machine tolerance)
    ## the sign is re-added
    if ( E %% 2 < .Machine$double.eps^0.5 )
      csim <- sgn*csim # restore sign
           #warning("E should be odd: ", E)
 
    ## 2: calculate total scoring S(i,c) and backtracing K(i,c)
    if ( verb>0 ) {
        cat(paste("Scoring matrix\t", time(), "\n",sep=""))
        cat(paste("parameters\t",paste("function:", S,
                                        "; scale:", E,
                                        "; max/min:", multi,sep=""),
                  "\n",sep=""))
    }
    ## TODO: handle Mn in scoring functions
    ## add official nuisance cluster
    SK<- calculateScore(seq=seqr, C=C, score=S, csim=csim,
                        M=M, Mn=Mn, multi=multi)

    ## 3: back-tracing to generate segments
    if ( verb>0 ) {
        cat(paste("Backtracing\t", time(), "\n",sep=""))
        cat(paste("parameters\t", paste("multib:",multib,sep=""), "\n",sep=""))
    }
    seg <- backtrace(S=SK$S, K=SK$K, multib=multib, nextmax=nextmax, verb=verb)

    ## 4: POST-PROCESSING
    ## remap: map back to original cluster names
    remap <- as.numeric(names(map))
    seg$segments[,1] <- remap[seg$segments[,1]]

    ## rm nuisance segments
    if ( rm.nui )
      seg$segments <- seg$segments[seg$segments[,1]!=0,,drop=FALSE]

    ## segmentation ID and sequence length
    seg$N <- N # sequence length
    seg$ids <- "segments" # default ID; required in plot functions
    
    ## add colors
    colors <- NULL
    ## inherit colors from cset
    if ( !is.null(cset) ) 
      if ( "colors" %in% names(cset) ) 
          colors <- cset$colors[[k]]
    ## input: plain sequence or cset without colors
    if ( is.null(colors) ) { 
        ## generate colors; use remap for names
        colors <- rep("#888888", length(remap)) # nuisance color!
        names(colors) <- sort(remap)
        colors[remap!=0] <- color_hue(length(remap[remap!=0])) # non-nuisance
    }
    seg$colors <- list(colors)
    names(seg$colors) <- seg$ids
    
    ## add matrices if requested!
    ## ... can be used for plotting or re-analysis
    if ( save.matrix ) {
        colnames(SK$S) <- colnames(SK$K) <- colnames(SK$S1) <- 0:(ncol(SK$S)-1)
        ## scoring matrices, S1, S, K
        seg$SK <- list(SK)
        ## add cluster similarity marix
        seg$csim <- list(csim) # not present in batch!
        names(seg$SK) <- names(seg$csim) <- seg$ids
    }    

    

    ## record run-time
    elapsed <- difftime(Sys.time(),stime,, units="secs")
    seg$elapsed <- elapsed

    ## add used parameters
    ## TODO: only parameters used for selected scoring
    parms <- data.frame(k=k,
                        S=S, E=E, M=M, Mn=Mn, a=a, nui=nui, multi=multi, 
                        nextmax=nextmax, multib=multib,
                        stringsAsFactors=FALSE)
    rownames(parms) <- seg$ids
    seg[["settings"]] <- parms


    ## assign S3 class
    class(seg) <- "segments"

    ## report timing
    if ( verb>0 ) {
        cat(paste("elapsed, sec\t", round(elapsed), "\n",sep=""))
        cat(paste("Done at  \t", time(), "\n",sep=""))
    }
    return(seg)
    
}

## TODO: move this to .cpp as well, then the whole algo is available in C++
#' Back-tracing step of the \code{segmenTier} algorithm.
#' 
#' back-tracing step: collect clustered segments from the scoring matrix
#' \code{S(i,c)} by back-tracing the position \code{j=k} which delivered
#' the maximal score at position \code{i}.
#' 
#' @param S matrix S, containing the local scores
#' @param K matrix K, containing the position k used for score maximization
#' @param multib if multiple k produce the maximal score, take either the
#' shortest k ("max") or the longest k ("min"); if \code{multib} is set to
#' "skip" the next unique k will be searched
#' @param nextmax proceed backwards while score is increasing before
#' opening a new segment
#' @param verb print messages
#' @export
backtrace <- function(S, K, multib, nextmax=FALSE, verb=TRUE) {

    segments <- NULL
    i <- nrow(S)
    warnings <- NULL
    
    while ( i>0 ) {
        
        ## FIND SEGMENT END
        ## Note, that this determines the segment's cluster assignment,
        ## unless multiple clusters deliver the maximal score
        
        ##  WHICH cluster(s) had max S at i?
        c <- which(S[i,]==max(S[i,]))
        
        ## search next max(S[i,c]) over i--
        if ( nextmax ) {
            while( i>0 & sum(S[i,c] <= S[i-1,c])==length(c) ) 
              i <- i - 1
            ##  which cluster had max S at i?
            c <- which(S[i,]==max(S[i,]))
        }

        ## FIND SEGMENT START
        ## WHICH k was used?

        ## get the k that was used for the maximum score S(i,c)
        k <- K[i,c]

        ## ... and handle multiple max scores from several clusters
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning:", length(c),
                       "clusters with maximal score; c:", paste(c,collapse=";"),
                       "with max_k:", paste(k,collapse=";"), ":",
                       ifelse(multib=="skip", paste("skipping breakpoint",i),
                              ifelse(multib=="max",
                                     "taking shorter segment",
                                     "taking longer segment")), "\n")
                                     
            warnings <- warn(w,warnings,verb=verb)
            if ( multib=="skip" ) {
                i <- i-1
                next 
            }
            ## multib: max is shortest k, min is longest k
            ## i.e. max delivers the shorter segment, min the longer segment 
            km <- get(multib,mode="function")(k)
            c <- c[which(k==km)] 
            k <- km
        }
        ## if we still can't find one cluster 
        if ( length(c)>1 ) {
            w <- paste(i,"back-trace warning STILL:", length(c),
                       "clusters with maximal score; c:",
                       paste(c,collapse=";"),"\n") 
            warnings <- warn(w,warnings,verb=verb)
            c <- c[1] #paste(c,collapse=";")
        }
        
        ## ignore k==i segments
        ## NOTE this should only happen for obsolete k<=i instead of k<i
        if ( i %in% k ) {
            w<- paste(i,"back-trace warning:",sum(k%in%i),"0 length segments\n")
            warnings <- warn(w,warnings,verb=verb)
            c <- c[k!=i]
            k <- k[k!=i]
        }
        if ( length(k)==0 ) {
            w <- paste(i,"back-trace warning: no k left\n")
            warnings <- warn(w,warnings,verb=verb)
            i <- i-1
            next
        }

        ## assign segment!
        segments <- rbind(segments, c(c,k,i))
        i <- k - 1
    }
    segments <- segments[order(as.numeric(segments[,2])),,drop=FALSE]
    if ( !is.null(segments) )
        colnames(segments) <- c("CL","start","end")
    list(segments=segments,warnings=warnings)
}

### UTILS

## get segment table from \code{\link{segmentClusters}}
## @param x result object returned by function \code{\link{segmentClusters}}
## @return Returns the segment table as a matrix
## @export
#result <- function(x) x$segments

#' Print method for segmentation result from \code{\link{segmentClusters}}.
#' @param x result object returned by function \code{\link{segmentClusters}}
#' @param ... further argument to \code{print.data.frame}
#' @export
print.segments <- function(x, ...) {
    cat(paste("\nSimilarity-based segmentation by dynamic programming:"))
    cat(paste("\nTotal length: ", x$N))
    cat(paste("\nSegments:\n"))
    print(x$segments, ...)
    cat(paste("\nParameters:\n"))
    print(x$settings)
    if ( "elapsed"%in%names(x) ) {
        tme <- x$elapsed
        cat(paste("\nRun time"))
        if ( length(tme)>1 ) {
            tme <- mean(tme)
            cat(paste(", average"))
        }
        cat(paste(": ", round(tme,3), "seconds\n"))
    }
}

    
### DATA SET DOC

#' Transcriptome time-series from budding yeast.
#' 
#' Transcriptome time-series data from a region encompassing
#' four genes and a regulatory upstream non-coding RNA in budding yeast.
#' The data set is described in more detail in the publication
#' Machne, Murray & Stadler (2017) <doi:10.1038/s41598-017-12401-8>.
#'
#' @name tsd
#' @docType data
NULL
