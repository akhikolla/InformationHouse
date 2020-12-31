

### DATA TRANSFORMATION UTILS
## get Discrete Fourier Transformation
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,seq_len(n),drop=FALSE]
    if ( n==1 )
        colnames(fft) <- "DC"
    else
        colnames(fft) <- c("DC",as.character(seq_len(n-1)))
    fft
}
## Fourier permutation
do.perm <- function(x, fft=NULL, perm, verb=0) {

    N <- ncol(x)
    if ( is.null(fft) ) fft <- get.fft(x)
    xam <- abs(fft)/N
    pvl <- matrix(0,nrow=nrow(fft), ncol=ncol(fft))
    dimnames(pvl) <- dimnames(fft)
    ## TODO: use apply and parallel!
    for ( i in seq_len(perm) ) {
        if ( verb>0 & i%%round(perm/10)==0 )
          cat(paste(round(i/perm,2)*100,"%, "))
        ## randomize columns and get Fourier
        rft <- get.fft(x[,sample(seq_len(ncol(x)))])
        ram <- abs(rft)/N
        pvl <- pvl + as.numeric(ram >= xam)
    }
    if ( verb ) cat("\n")
    Re(pvl/perm)
}

#' \code{asinh} data transformation
#'
#' The asinh transformation, (\code{ash(x) = log(x + sqrt(x^2+1))}), is
#' an alternative to log transformation that has less (compressing) effects
#' on the extreme values (low and high values), and naturally handles
#' negative numbers and 0. Also see \code{\link{log_1}}.
#' @param x a numeric vector
#' @export
ash <- function(x) log(x+sqrt(x^2+1))

#' log transformation handling zeros by adding 1
#'
#' A conventional approach to handle 0 in log transformation is to simply
#' add 1 to all data, \code{log_1(x) = log(x+1)}. Also see \code{\link{ash}}.
#' @param x a numeric vector
#' @export
log_1 <- function(x) log(x+1)

## moving average
ma <- function(x, n=5, circular=FALSE) {
    stats::filter(x,rep(1/n,n), sides=2, circular=circular)
}

# calculate 95% confidence intervals for the given
# data vector using a t-distribution
ci95 <- function(data,na.rm=FALSE) {
    if ( na.rm ) data <- data[!is.na(data)]
    n <- length(data)
    if ( n<2 ) return(NA)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}

## cluster/segment colors; function derived from scale_colour_hue in ggplot2
## TODO: describe better
color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}


#' Process a time-series for clustering and segmentation.
#' 
#' Prepares a time-series (time points in columns) for subsequent
#' clustering, and performs requested data transformations, including
#' a Discrete Fourier Transform (DFT) of the time-series, as direct
#' input for the clustering wrapper
#' \code{\link{clusterTimeseries}}. When used for segmentation
#' the row order reflects the order of the data points along which
#' segmentation will occur. The function can also be used as a
#' stand-alone function equipped especially for analysis of
#' oscillatory time-series, including calculation of phases and
#' p-values for all DFT components, and can also be used for
#' Fourier Analysis and subsequent clustering without segmentation.
#' 
#' @details This function exemplifies the processing of an oscillatory
#' transcriptome time-series data as used in the establishment of this
#' algorithm and the demo \code{segment_data}. As suggested by Machne & Murray
#' (PLoS ONE 2012) and Lehmann et al. (BMC Bioinformatics 2014) a Discrete
#' Fourier Transform of time-series data allows to cluster time-series by
#' their change pattern. 
#'
#' Note that NA values are here interpreted as 0. Please take care of NA
#' values yourself, if you do not want this behavior.
#'
#' Rows consisting only of 0 (or NA) values, or with a total signal
#' (sum over all time points) below the value passed in argument
#' \code{low.thresh}, are detected, result in NA values in the
#' transformed data, and will be assigned to the
#' "nuisance" cluster in \code{\link{clusterTimeseries}}.
#'
#' Discrete Fourier Transform (DFT): if requested (option
#' \code{use.fft=TRUE}), a DFT will be applied using base R's
#' \code{\link[stats:fft]{mvfft}} function and reporting all or only
#' requested (option \code{dft.range}) DFT components, where the
#' first, or DC ("direct current") component, equals the total signal
#' (sum over all points) and other components are numbered 1:n,
#' reflecting the number of full cycles in the time-series. Values are
#' reported as complex numbers, from which both amplitude and phase
#' can be calculated.  All returned DFT components will be used by
#' \code{\link{clusterTimeseries}}.
#'
#' Additional Transformations: data can be transformed prior to DFT
#' (options \code{trafo}, \code{smooth.time}, \code{smooth.space}), or
#' after DFT (options \code{use.snr} and \code{dc.trafo}). It is
#' recommended to use the amplitude scaling (a signal-to-noise ratio
#' transformation, see option documentation).  The separate
#' transformation of the DC component allows to de-emphasize the total
#' signal in subsequent clustering & segmentation.  Additionally, but
#' not tested in the context of segmentation, a Box-Cox transformation
#' of the DFT can be performed (option \code{lambda}).  This
#' transformation proofed useful in DFT-based clustering with the
#' model-based clustering algorithm in package \pkg{flowClust}, and is
#' available here for further tests with k-means clustering.
#' 
#' Phase, Amplitude and Permutation Analysis: this time-series
#' processing and subsequent clustering can also be used without
#' segmentation, eg. for conventional microarray data or RNA-seq data
#' already mapped to genes. The option \code{perm} allows to perform a
#' permutation test (\code{perm} times) and adds a matrix of empirical
#' p-values for all DFT components to the results object, ie. the
#' fraction of \code{perm} where amplitude was higher then the
#' amplitude of the randomized time-series.  Phases and amplitudes can
#' be derived from the complex numbers in matrix "dft" of the result
#' object.
#' @param ts a time-series as a matrix, where columns are the
#'     time points and rows are ordered measurements, e.g., genomic
#'     positions for transcriptome data
#' @param na2zero interpret NA values as 0
#' @param trafo prior data transformation, pass any function name,
#'     e.g., "log", or the package functions "ash" (asinh:
#'     \code{ash(x) = log(x + sqrt(x^2+1))}) or "log_1"
#'     (\code{log(ts+1)})
#' @param low.thresh use this threshold to cut-off data, which will be
#'     added to the absent/nuisance cluster later
#' @param perm number of permutations of the data set, to obtain
#'     p-values for the oscillation
#' @param use.fft use the Discrete Fourier Transform of the data
#' @param dft.range a vector of integers, giving the components of the
#'     Discrete Fourier Transform to be used where 1 is the first
#'     component (DC) corresponding to the total signal (sum over all
#'     time points), and 2:n are the higher components corresponding
#'     to 2:n full cycles in the data
#' @param use.snr use a scaled amplitude, where each component of the
#'     Discrete Fourier Transform is divided by the mean of all other
#'     components (without the first or DC component), a normalization
#'     that can be interpreted to reflect a signal-to-noise ratio
#'     (SNR)
#' @param lambda parameter lambda for Box-Cox transformation of DFT
#'     amplitudes (experimental; not tested)
#' @param dc.trafo data transformation for the first (DC) component of
#'     the DFT, pass any function name, e.g., "log", or the package
#'     functions "ash" (asinh: \code{ash(x) = log(x + sqrt(x^2+1))})
#'     or "log_1" (\code{log(x+1)}).
#' @param smooth.space integer, if set a moving average is calculated
#'     for each time-point between adjacent data points using stats
#'     package's \code{\link[stats:smooth]{smooth}} with option
#'     \code{span=smooth.space}
#' @param smooth.time integer, if set the time-series will be smoothed
#'     using stats package's \code{\link[stats:filter]{filter}} to
#'     calculate a moving average with span \code{smooth.time} and
#'     \code{\link[stats:smoothEnds]{smoothEnds}} to extrapolate
#'     smoothed first and last time-points (again using span
#'     \code{smooth.time})
#' @param circular.time logical value indicating whether time can be
#'     treated as circular in smoothing via option \code{smooth.time}
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @return Returns a list of class "timeseries" which comprises of
#'     the transformed time-series and additional information, such as
#'     the total signal, and positions of rows with only NA/0
#'     values. Note that NA values are interpreted as 0.
#' @references Machne & Murray (2012)
#'     <doi:10.1371/journal.pone.0037906>, and Lehmann et al. (2013)
#'     <doi:10.1186/1471-2105-14-133>
#' @examples
#' data(primseg436)
#' ## The input data is a matrix with time points in columns
#' ## and a 1D order, here 7624 genome positions, is reflected in rows,
#' ## if the time-series should be segmented.
#' nrow(tsd)
#' ## Time-series processing prepares the data for clustering,
#' ## the example data is periodic, and we will cluster its Discrete Fourier
#' ## Transform (DFT) rather then the original data. Specifically we will
#' ## only use components 1 to 7 of the DFT (dft.range) and also apply
#' ## a signal/noise ratio normalization, where each component is
#' ## divided by the mean of all other components. To de-emphasize
#' ## total levels the first component (DC for "direct current") of the
#' ## DFT will be separately arcsinh transformed. This peculiar combination
#' ## proofed best for our data:
#' tset <- processTimeseries(ts=tsd, na2zero=TRUE, use.fft=TRUE,
#'                           dft.range=1:7, dc.trafo="ash", use.snr=TRUE)
#' ## a plot method exists for the returned time-series class:
#' par(mfcol=c(2,1))
#' plot(tset)
#'@export
processTimeseries <- function(ts, na2zero=FALSE, trafo="raw", 
                              use.fft=FALSE, dc.trafo="raw", dft.range,
                              perm=0, use.snr=FALSE, lambda=1,
                              low.thresh=-Inf, 
                              smooth.space=1,
                              smooth.time=1, circular.time=FALSE,
                              verb=0) {

    if ( typeof(ts)=="list" )
        ts <- as.matrix(ts) # smoothEnds causes problems for data.frames!?
     
    
    tsd <-ts
    ## NOTE: replace NA by 0
    ## TODO: make this behavior optional?
    if ( na2zero ) {
        if ( any(is.na(tsd)) )
            warning("Setting ", sum(is.na(tsd)) , " NA values to 0.")
        tsd[is.na(tsd)] <- 0 # set NA to zero (will become nuisance cluster)
    }
    if ( any(is.na(tsd)) )
        warning(sum(is.na(tsd)) , " NA values in input data detected.")
    
    ## detect rows only consisting of 0, these will not be processed
    ## and later assigned to a nuisance cluster
    zs <- apply(tsd==0,1,sum)==ncol(tsd) # remember all zeros

    ## smooth time-points between adjacent positions
    ## currently not used, doesn't help to avoid fragmentation!
    if ( smooth.space>1 ) {
        tsm <- apply(tsd[!zs,], 2 ,ma, smooth.space,FALSE)
        tsd[!zs,] <- tsm
    }
    ## smooth time-series
    ## currently used only in clustering final segment time series
    ## NOTE/TODO: smooth.time must be ODD for smoothEnds
    if ( smooth.time>1 ) {
        tsm <- t(apply(tsd[!zs,], 1, ma, n=smooth.time, circular=circular.time))
        ends <- stats::smoothEnds(tsd[!zs,], k=smooth.time)
        tsd[!zs,] <- tsm
        tsd[!zs,c(1,ncol(tsd))] <- ends[,c(1,ncol(tsd))]
    }

        
    ## transform raw data?
    ## NOTE that DFT and SNR below (use.fft) are an alternative
    ## data normalization procedure
    ## default: identity
    if ( trafo!="raw" )
        tsd <- get(trafo, mode="function")(tsd) # ash, log_1, etc
    
    ## get DFT
    if ( missing(dft.range) )
        dft.range <- NULL
    fft <- pvl <- NULL
    if ( use.fft ) {

        ## get DFT
        tmp <- get.fft(tsd[!zs,])
        fft <- matrix(NA, ncol=ncol(tmp), nrow=nrow(tsd))
        rownames(fft) <- rownames(tsd)
        colnames(fft) <- colnames(tmp)
        fft[!zs,] <- tmp
        

        ## do DFT on permuted time-series to obtain p-values
        ## TODO: include amplitude and DC scaling in permutation?
        if ( perm>0 ) { 
            tmp <- do.perm(tsd[!zs,],fft=fft[!zs,], perm, verb=verb)
            pvl <- matrix(NA, ncol=ncol(tmp), nrow=nrow(tsd))
            colnames(pvl) <- colnames(tmp)
            pvl[!zs,] <- tmp
        }
        
        ## amplitude-scaling (~SNR), see Machne&Murray 2012
        if ( use.snr ) {
            amp <- abs(fft)
          snr <- fft
          for ( a in 2:ncol(fft) )
            snr[,a] <- fft[,a]/apply(amp[,-c(1,a)],1,mean)
          fft <- snr
        }

        ## experimental: amplitude Box-Cox transformation
        ## box-cox trafo for negative values (Bickel and Doksum 1981)
        ## as used in flowClust
        bc <- function(x,lambda) (sign(x)*abs(x)^lambda-1)/lambda
        ## amplitude box-cox trafo for complex polar coordinates 
        bcdft <- function(x, lambda) {
            if ( class(x)=="matrix" )
                return(apply(x,2, bcdft, lambda))
            ## Box-Cox transform amplitude
            y <- bc(abs(x), lambda)
            ## amplitude scaling factor
            sf <- (y-min(y,na.rm=TRUE))/abs(x)
            x*sf
        }
        if ( lambda!=1 ) {
            fft <- bcdft(fft, lambda=lambda)
        }
        
        ## PREPARE DATA FOR CLUSTERING
        ## get low expression filter!
        tot <- Re(fft[,1]) # NOTE: DC component = rowSums(tsd)
        low <- tot < low.thresh

        ## DC scaling
        if ( dc.trafo!="raw" )
            fft[,1] <- get(dc.trafo,mode="function")(fft[,1]) # ash, log_1, etc

        ## filter selected components
        if ( is.null(dft.range) ) # allows passing NULL to use auto
            dft.range <- seq_len(ncol(fft))
        dat <- fft[,dft.range,drop=FALSE]
         
        ## get Real and Imaginary pars
        re <- Re(dat)
        colnames(re) <- paste("Re_",colnames(re),sep="")
        im <- Im(dat)
        colnames(im) <- paste("Im_",colnames(im),sep="")
        #if ( 1 %in% dft.range ) # rm 0 DC component from Im
        #    im <- im[,-1]
        ## filter 0 imaginary components: DC and Nyquist!
        im <- im[,apply(im,2,function(x) any(x!=0,na.rm=TRUE))]
        dat <- cbind(re,im)

    }else {

        dat <- tsd
        dat[zs,] <- NA # set zero-vals to NA

        ## get low expression filter
        tot <- rowSums(dat,na.rm=TRUE)
        low <- rep(FALSE, nrow(dat))
        low <- tot < low.thresh
    }

    ## NA can come from init of fft matrix
    ## complete time series are NA 
    na.rows <- rowSums(is.na(dat))==ncol(dat)

    ## only some fields are NA
    na.fields <- is.na(rowSums(dat,na.rm=FALSE)) 
#    chk <- sum(na.fields & !na.rows)
#    if ( chk>0 )
#        warn(chk, " rows have individual NA fields and are removed;",
#             "please manually set those if they are to be retained")
    
    ## check 0 variance
    ##no.var <- apply(dat, 1, var)==0

    ## remove data rows: NA or low, or 0 variance
    rm.vals <- (na.rows|na.fields) | (low)# | no.var) 

    settings <- list(trafo=trafo, 
                     use.fft=use.fft,
                     dc.trafo=dc.trafo,
                     dft.range=dft.range,
                     lambda=lambda,
                     perm=perm,
                     use.snr=use.snr,
                     low.thresh=low.thresh, 
                     smooth.space=smooth.space,
                     smooth.time=smooth.time)

    ## generate processing ID - this will be inherited to clusters
    ## and from there to segment ID and type
    processing <- paste("T:",trafo,sep="")
    if ( use.fft )
      processing <- paste(processing,"_",
                          paste("D:dft",paste(range(dft.range),collapse="-"),
                                sep=""),".",
                          paste("dc",dc.trafo,sep=""),".",
                          ifelse(use.snr,"snr","raw"),
                          sep="")

    ## time-series data set for clustering in clusterTimeseries
    tset <- list(dat=dat, ts=tsd, dft=fft, pvalues=pvl, tot=tot,
                 zero.vals=zs, rm.vals=rm.vals, low.vals=low,
                 settings=settings, id=processing)
    class(tset) <- "timeseries"
    
    ## silent return
    tmp <- tset
}

## TODO: adapt to be used in segmentation as well, is fcls@mu
## equal/similar to kmeans' 'centers'? Are Ccc and Pci calculated correctly?

#' Cluster a processed time-series with
#' \code{\link[flowClust:flowClust]{flowClust}} &
#' \code{\link[flowMerge:merge]{flowMerge}}.
#' 
#' A wrapper for \code{\link[flowClust:flowClust]{flowClust}}, clustering
#' a time-series object \code{tset} provided by \code{\link{processTimeseries}},
#' where specifically the DFT of a time-series and requested data
#' transformation were calculated. This is intended to work in the same way
#' as \code{\link{clusterTimeseries}} but was so far only tested for
#' clustering of the final segment time-series, as previously applied
#' to microarray data from yeast by Machne & Murray (2012)
#' <doi:10.1371/journal.pone.0037906> and from cyanobacteria by Lehmann
#' et al. (2013) <doi:10.1186/1471-2105-14-133>.
#' It could in principle also be used for segmentation, but that has not
#' been extensively tested. \code{\link[flowClust:flowClust]{flowClust}}
#' implements a model-based clustering approach and is much slower then
#' \code{\link[stats:kmeans]{kmeans}} used in \code{\link{clusterTimeseries}}. 
#' Please see option \code{ncpu} on how to use parallel mode, which
#' does not work on some installations. However, model-based clustering has
#' the advantage of an intrinsic measure (\code{BIC}) to decide on the optimal
#' cluster numbers. Additionally, the clusters can be "merged" to fewer
#' clusters at constant \code{BIC} using
#' \code{\link[flowMerge:merge]{flowMerge}}.
#' @param tset processed time-series as provided by
#' \code{\link{processTimeseries}}
#' @param ncpu number of cores available for parallel mode of
#' \pkg{flowClust}. NOTE: parallel mode of
#' \code{\link[flowClust:flowClust]{flowClust}} is often non-functional.
#' Alternatively, you can set \code{options(mc.cores=ncpu)} directly.
#' @param K the requested cluster numbers (vector of integers)
#' @param merge logical indicating whether cluster merging with
#'  \code{\link[flowMerge:merge]{flowMerge}} should be attempted
#' @param selected a pre-selected cluster number  which is then
#' used as a start clustering for  \code{\link[flowMerge:merge]{flowMerge}}
#' (if option \code{merge==TRUE})
#' @param B maximal number of EM iterations 
#' @param tol tolerance for EM convergence
#' @param lambda initial Box-Cox trafo
#' @param nu degrees of freedom used for the t distribution, Inf for
#' pure Gaussian
#' @param nu.est 0: no, 1: non-specific, 2: cluster-specific estimation of nu
#' @param trans 0: no, 1: non-specific, 2: cluster-specific estim. of lambda
#' @param ... further parameters for
#' \code{\link[flowClust:flowClust]{flowClust}}
#' @references Machne & Murray (2012)
#'     <doi:10.1371/journal.pone.0037906>
#' @export
flowclusterTimeseries <- function(tset, ncpu=1, K=10, selected, merge=FALSE,
                                  B=500, tol=1e-5, lambda=1,
                                  nu=4, nu.est=0, trans=1, ...) {

    if ( !requireNamespace("flowClust", quietly = TRUE) )
        stop("`flowclusterTimeseries' requires the bioconductor package ",
             "`flowClust'")
    if ( merge & !requireNamespace("flowMerge", quietly = TRUE) )
      stop("option `merge' requires the bioconductor package `flowMerge'")
    
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    clsDat <- dat[!rm.vals,]

    ##if ( ncpu>1 ) {
    ## TODO: when is which required?
        oldcpu <- unlist(options("cores"))
        oldcpu2 <- unlist(options("mc.cores"))
        options(cores=ncpu)
        options(mc.cores=ncpu)
    ##}

    ## NOTE: passing mc.cores=ncpu to flowClust (via ...)
    ## works on some installations!
    fcls <- flowClust::flowClust(clsDat, K=K, B=B, tol=tol, lambda=lambda,
                                 nu=nu, nu.est=nu.est, trans=trans,...)

    ## collect clusterings
    cluster.matrix <- matrix(0, nrow=nrow(dat), ncol=length(K))
    colnames(cluster.matrix) <- as.character(K)
    rownames(cluster.matrix) <- rownames(dat)
    bic <- rep(NA, length(K))
    names(bic) <- as.character(K)
    icl <- bic
    for ( i in seq_along(fcls) ) {
      if ( length(fcls) > 1 ) fc <- fcls[[i]]
      else fc <- fcls
      cl.num <- as.character(fc@K)
      cluster <- flowClust::Map(fc,rm.outliers=FALSE)
      cluster.matrix[!rm.vals, cl.num] <- cluster
      bic[cl.num] <- fc@BIC
      icl[cl.num] <- fc@ICL
    }
    ## max BIC and ICL
    max.bic <- max(bic, na.rm=TRUE)
    max.clb <- K[which(bic==max.bic)]
    max.icl <- max(icl, na.rm=TRUE)
    max.cli <- K[which(icl==max.icl)]
    ## best K selection
    ## use K with max BIC, unless specified in argument
    if ( missing(selected) )
      selected <- max.clb
    if ( selected==0 )
      selected <- max.clb

    ## MERGE CLUSTERS, starting from best BIC by flowMerge
    mrg.orig <- mrg.cl <- mrg.id <-  mrg.nm <- obj <- NULL
    if ( merge ) {
        best <- which(K==selected)
        if ( length(fcls) > 1 ) fc <- fcls[[best]]
        else fc <- fcls
        mcls <- rep(0, nrow(dat))
        mrg.id <- mrg.nm <- mrg.cl <- "NA"
        obj <- try(flowMerge::flowObj(fc, flowCore::flowFrame(clsDat)))
        if ( class(obj)!="try-error" ) {
            mrg <- try(flowMerge::merge(obj))
            if ( class(mrg)!="try-error" ) {
                mrg.cl <- flowMerge::fitPiecewiseLinreg(mrg)
                obj <- mrg[[mrg.cl]]
                mcls[!rm.vals] <- flowClust::Map(obj, rm.outliers=FALSE)
                mrg.orig <- K[best] # source K
                mrg.id <- paste0(K[best],"m",mrg.cl) # merged K
                mrg.nm <- paste0("K:",K[best],"m",mrg.cl) # final column name
            }
        }
        cluster.matrix <- cbind(cluster.matrix, mcls)
        colnames(cluster.matrix)[ncol(cluster.matrix)] <- mrg.id
        all <- append(fcls,obj)
    } else {
        all <- fcls
    }
    colnames(cluster.matrix) <- paste0("K:",colnames(cluster.matrix))

    ## collect centers, Pci and Ccc correlation matrices (see clusterTimeseries)
    ## TODO: TEST FOR SEGMENTATION (currently only used for final
    ## segment time series)
    ## -> is `mu' really the same as centers and are Ccc and Pci
    ## correct? 
    centers <- Pci <- Ccc <- rep(list(NA),length(all))
    for ( i in seq_along(all) ) {
        if ( length(all) > 1 ) fc <- all[[i]]
        else fc <- all

        ## get cluster centers!
        x <- fc@mu
        rownames(x) <- seq_len(nrow(x))
        colnames(x) <- colnames(clsDat)
        
        centers[[i]] <- x
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(x))

        Ccc[[i]] <- cr
        
        ## P(c,i) - position X cluster correlation
        P <- matrix(NA,nrow=nrow(dat),ncol=nrow(fc@mu))
        P[!rm.vals,] <- clusterCor_c(clsDat, fc@mu)

        Pci[[i]] <- P
    }
    names(centers) <- names(Pci) <- names(Ccc) <- colnames(cluster.matrix)

    
    ## clustering data set for use in segmentCluster.batch 
    fcset <- list(clusters=cluster.matrix,
                  N=sum(!rm.vals), # number of clustered data
                  M=ncol(dat), # dimension of clustered data
                  centers=centers, Pci=Pci, Ccc=Ccc,
                  K=K, usedk=K, selected=selected, warn=NULL,
                  bic=bic, icl=icl,
                  ids=colnames(cluster.matrix),
                  tsid=rep(tset$id,ncol(cluster.matrix)),
                  flowClust=fcls, flowMerge=obj, # flowClust/flowMerge results
                  max.clb=max.clb, max.cli=max.cli,
                  merged.origK=mrg.orig, merged.K=mrg.cl, merged=mrg.nm)
    class(fcset) <- "clustering" 

    ## add cluster colors
    fcset <- colorClusters(fcset)

    ##if ( ncpu>1 )
    options(cores=oldcpu)
    options(mc.cores=oldcpu2)
    
    ## silent return
    tmp <- fcset
}

#' Experimental: AIC/BIC for kmeans
#'
#' This function is supposed to provide a log-likelihood method for
#' \code{\link[stats:kmeans]{kmeans}} results, after Neal Fultz at
#' \url{https://stackoverflow.com/a/33202188} and also featured in the
#' \href{https://rdrr.io/github/nfultz/stackoverflow/src/R/logLik_kmeans.R}{stackoverflow package}. Note, that the blogged version on Jan 30, 2019 adds
#' a minus and a division by 2 compared to a linked git version.
#' This idea has not been reviewed, and this function has not been
#' tested extensively; feel free to do so and contribute your results.
#' 
#' This is an attempt to reproduce the \code{BIC} measure
#' in model-based clustering to decide on an optimal number of clusters.
#' This function will be used for \code{\link[stats:kmeans]{kmeans}}
#' results objects when passed to \code{\link[stats:BIC]{BIC}} and
#' \code{\link[stats:AIC]{AIC}} functions from the \pkg{stats} package in
#' base R, and BIC and AIC are calculated this way in
#' \code{\link{segmentClusters}}. It is however not used anywhere at the
#' moment.
#' @param object a \code{\link[stats:kmeans]{kmeans}} result object
#' @param ... unused
#' @export
logLik.kmeans <- function(object, ...)
    structure(-object$tot.withinss/2,
              df = nrow(object$centers)*ncol(object$centers),
              nobs = length(object$cluster),
              class = 'logLik')


#' Cluster a processed time-series with k-means.
#' 
#' Performs \code{\link[stats:kmeans]{kmeans}} clustering of a
#' time-series object \code{tset} provided by
#' \code{\link{processTimeseries}}, and calculates cluster-cluster
#' and cluster-position similarity matrices as required for
#' \code{\link{segmentClusters}}.
#'
#' @details This function performs one or more time-series clustering(s)
#' using \code{\link[stats:kmeans]{kmeans}}, and the output of
#' \code{\link{processTimeseries}} as input. It further calculates
#' cluster centers, cluster-cluster and cluster-position similarity
#' matrices (Pearson correlation) that will be used by the main function
#' of this package, \code{\link{segmentClusters}}, to split the cluster
#' association sequence into segments, and assigns each segment to
#' the "winning" input cluster.
#'
#' The argument \code{K} is an integer vector that sets the requested
#' cluster numbers (argument \code{centers} in
#' \code{\link[stats:kmeans]{kmeans}}). However, to avoid errors in batch
#' use, a smaller \code{K} is chosen, if the data contains less then
#' \code{K} distinct values.
#' 
#' Nuisance Cluster:
#' values that were removed during time-series processing, such as
#' rows that only contain 0 or NA values, will be assigned to
#' the "nuisance cluster" with cluster label "0". Additionally, a minimal
#' correlation to any cluster center can be specified, argument
#' \code{nui.thresh}, and positions without any correlation higher
#' then this, will also be assigned to the "nuisance" cluster.
#' Resulting "nuisance segments" will not be shown in the results.
#'
#' Cluster Sorting and Coloring:
#' additionally the cluster labels in the result object will be sorted by
#' cluster-cluster similarity (see \code{\link{sortClusters}}) and cluster
#' colors assigned (see \code{\link{colorClusters}}) for convenient data
#' inspection with the plot methods available for each data processing
#' step (see examples).
#' 
#' Note that the function, in conjunction with
#' \code{\link{processTimeseries}}, can also be used as a stand-alone
#' tool for time-series clusterings, specifically implementing the
#' strategy of clustering the Discrete Fourier Transform of periodic
#' time-series developed by Machne & Murray (2012)
#' <doi:10.1371/journal.pone.0037906>, and further analyzed in Lehmann
#' et al. (2013) <doi:10.1186/1471-2105-14-133>, such as transcriptome
#' data from circadian or yeast respiratory oscillation systems.
#' @param tset a "timeseries" object returned by
#'     \code{\link{processTimeseries}}
#' @param K the number of clusters to be calculated, ie. the argument
#'     \code{centers} of \code{\link[stats:kmeans]{kmeans}}, but here
#'     multiple clusterings can be calculated, ie. \code{K} can be an
#'     integer vector. Note that a smaller cluster number is automatically
#'     chosen, if the data doesn't have more then K different values.
#' @param iter.max the maximum number of iterations allowed in
#'     \code{\link[stats:kmeans]{kmeans}}
#' @param nstart number of randomized initializations of
#'     \code{\link[stats:kmeans]{kmeans}}: "how many random sets should
#'     be chosen?"
#' @param nui.thresh threshold correlation of a data point to a
#'     cluster center; if below the data point will be added to
#'     nuisance cluster 0
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @return Returns a list of class "clustering" comprising of a matrix
#'     of clusterings, lists of cluster centers, cluster-cluster and
#'     cluster-position similarity matrices (Pearson correlation) used
#'     by \code{\link{segmentClusters}}, and additional information
#'     such as a cluster sorting by similarity and cluster colors that
#'     allow to track clusters in plots. A plot method exists that
#'     allows to plot clusters aligned to "timeseries" and "segment"
#'     plots.
#' @references Machne & Murray (2012)
#'     <doi:10.1371/journal.pone.0037906>, and Lehmann et al. (2013)
#'     <doi:10.1186/1471-2105-14-133>
#' @examples
#' data(primseg436)
#' ## Discrete Fourier Transform of the time-series, 
#' ## see ?processTimeseries for details
#' tset <- processTimeseries(ts=tsd, na2zero=TRUE, use.fft=TRUE,
#'                           dft.range=1:7,  dc.trafo="ash", use.snr=TRUE)
#' ## ... and cluster the transformed time-series
#' cset <- clusterTimeseries(tset)
#' ## plot methods for both returned objects allow aligned plots
#' par(mfcol=c(3,1))
#' plot(tset)
#' plot(cset)
#'@export
clusterTimeseries <- function(tset, K=16, iter.max=100000, nstart=100,
                              nui.thresh=-Inf, verb=1) {


    ## TODO:
    ## call recursively if multiple tsets are available
    ## and pre-pend tset names 

    ## get time series data
    id <- tset$id
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    N <- nrow(dat)

    ## enough distinct values?
    ## TODO: issue segment based on low-filter
    ## OOR: postprocessing - extend segments into low levels?
    warn <- NULL
    if ( sum(!rm.vals)<10 ) 
        warn <- "not enough data"
    else if ( sum(!duplicated(dat[!rm.vals,]))<2 ) 
        warn <- "not enough data diversity"
    if ( !is.null(warn) ) {
        warning(warn)
        return(NULL)
    }
    
    ## CLUSTERING
    ## stored data
    clusters <- matrix(NA, nrow=nrow(dat), ncol=length(K))
    rownames(clusters) <- rownames(dat)
    centers <- Pci <- Ccc <- rep(list(NA), length(K))

    ## BIC/AIC for kmeans
    bic <- rep(NA, length(K))
    names(bic) <- as.character(K)
    aic <- bic
    
    if ( verb>0 ) {
        cat(paste("Timeseries N\t",N,"\n",sep=""))
        cat(paste("Used datapoints\t",sum(!rm.vals),"\n",sep=""))
    }
    
    usedk <- K
    for ( k in seq_along(K) ) {
        
        ## get cluster number K
        ## NOTE: a smaller cluster number is automatically chosen,
        ## if the data doesn't have more then K different values!
        Kused <- min(c(K[k],sum(!duplicated(dat[!rm.vals,]))))
        
        if ( verb>0 )
            cat(paste("Clusters K\t", Kused, "\n",sep=""))
        
        ## cluster
        km <- stats::kmeans(dat[!rm.vals,], Kused, iter.max=iter.max,
                            nstart=nstart, algorithm="Hartigan-Wong")
        ## use alternative algo if this error occurred
        if (km$ifault==4) {
            km <- stats::kmeans(dat[!rm.vals,], Kused,
                                iter.max=iter.max,nstart=nstart,
                                algorithm="MacQueen")
            warn <- paste("quick-transfer error in kmeans algorithm",
                          "Hartigan-Wong, taking MacQueen")
            warning(warn)
        }
        
        ## prepare cluster sequence
        seq <- rep(0, N) ## init. to nuisance cluster 0
        seq[!rm.vals] <- km$cluster
        
        ## store which K was used, the clustering and cluster centers
        usedk[k] <- Kused
        clusters[,k] <- seq
        centers[[k]] <- km$centers
        
        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(km$centers))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster correlation
        P <- matrix(NA,nrow=N,ncol=Kused)
        P[!rm.vals,] <- clusterCor_c(dat[!rm.vals,], km$centers)

        Pci[[k]] <- P

        ## calculate BIC/AIC
        ## NOTE: this uses logLik.kmeans defined herein!
        bic[k] <- stats::BIC(km)
        aic[k] <- stats::AIC(km)
    }

    ## re-assign by correlation threshold
    ## NOTE: this only affects scoring function ccor
    for ( k in seq_len(ncol(clusters)) ) {
        cls <- clusters[,k]
        for ( p in seq_len(nrow(Pci[[k]])) )
            if ( !any(Pci[[k]][p,] > nui.thresh, na.rm=TRUE) )
                cls[p] <- 0
        clusters[,k] <- cls
    }           
    
    ## count duplicate K
    if ( any(duplicated(K)) ) {
        sel <- paste(K,".1",sep="")
        cnt <- 2
        while( sum(duplicated(sel)) ) {
            sel[duplicated(sel)] <- sub("\\..*",paste0(".",cnt),
                                        sel[duplicated(sel)])
            cnt <- cnt+1
        }
        K <- sub("\\.1$","",sel)
    }
    ## name all results by K, will be used!
    colnames(clusters) <- names(centers) <-
        names(Pci) <- names(Ccc) <- paste0("K:",K) #paste(id,"_K:",K,sep="")

    ## max BIC and ICL
    max.bic <- max(bic, na.rm=TRUE)
    max.clb <- K[which(bic==max.bic)[1]]
    max.aic <- max(aic, na.rm=TRUE)
    max.cla <- K[which(aic==max.aic)[1]]
    ## best K selection
    ## use K with max BIC
    selected <- max.clb
    
    ## clustering data set for use in segmentCluster.batch 
    cset <- list(clusters=clusters, 
                 N=sum(!rm.vals), # number of clustered data
                 M=ncol(dat), # dimension of clustered data
                 centers=centers, Pci=Pci, Ccc=Ccc,
                 K=K, usedk=usedk, selected=selected,
                 bic=bic, aic=aic, 
                 max.clb=max.clb, max.cla=max.cla,
                 warn=warn, ids=colnames(clusters),
                 tsid=rep(id,ncol(clusters)))
    class(cset) <- "clustering"

    ## add cluster colors
    cset <- colorClusters(cset)
    
    ## silent return
    tmp <- cset
}

#' Assign colors to clusters.
#' 
#' Takes a clustering set as returned by \code{\link{clusterTimeseries}} and
#' assigns colors to each cluster in each clustering along
#' the "hue" color wheel, as in \code{scale_colour_hue} in \code{ggplot2}.
#' If \code{cset} contains a sorting, this sorting will be used to assign
#' colors along the color wheel, otherwise a sorting will be calculated first,
#' using \code{\link{sortClusters}}.
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#' @param colf a function that generates \code{n} colors
#' @param ... arguments to color function \code{colf}
#' @return Returns the input "clustering" object with a list of vectors
#' ("colors"), each providing a named vector of colors for each cluster.
#'@export
colorClusters <- function(cset, colf, ...) {

    ## each column in the clustering matrix is one clustering
    cset$colors <- rep(list(NA), ncol(cset$clusters))

    ## sort if no sorting is present
    ## this requires the cluster-cluster similarity matrix Ccc
    if ( !"sorting" %in% names(cset) )
        cset <- sortClusters(cset, sort=TRUE)

    ## generate colors; use gray for nuisance
    for ( k in seq_len(ncol(cset$clusters)) ) {
        if ( missing(colf) )
            colf <- color_hue # internal function
        srt <- cset$sorting[[k]]
        srt <- srt[srt!="0"]
        cols <- c("#888888", colf(length(srt)),...)
        names(cols) <- c("0", srt)
        cset$colors[[k]] <- cols
    }
    names(cset$colors) <- colnames(cset$clusters)
    cset
}

#' Sort clusters by similarity.
#' 
#' Takes a "clustering" object as returned by
#' \code{\link{clusterTimeseries}} and uses the cluster-cluster
#' similarity matrix, item \code{Ccc}, to sort clusters by their
#' similarity, starting with the cluster labeled `1'; the next
#' cluster is the first cluster (lowest cluster label) with the
#' highest similarity to cluster `1', and proceeding from there. The
#' final sorting is added as item "sorting" to the \code{cset} object
#' and returned.  This sorting is subsequently used to select cluster
#' colors and in the plot method. This simply allows for more
#' informative plots of the clustering underlying a segmentation but
#' has no consequence on segmentation itself.
#' @param cset a clustering set as returned by
#'     \code{\link{clusterTimeseries}}
#' @param sort if set to FALSE the clusters will be sorted merely
#'     numerically
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @return Returns the input "clustering" object with a list of
#'     vectors (named "sorting"), each providing a similarity-based
#'     sorting of cluster labels.
#'@export
sortClusters <- function(cset, sort=TRUE, verb=0) {

    ## each column in the clustering matrix is one clustering
    cset$sorting <- rep(list(NA), ncol(cset$clusters))

    ## merely generate numerical sorting if sort is FALS
    if ( !sort ) {
        sorting <- NULL
        for ( k in seq_len(ncol(cset$clusters)) ) 
            sorting[[k]] <- sort(unique(cset$clusters[,k]))
        sorting <- lapply(sorting, function(x) x[x!=0])
        names(sorting) <- colnames(cset$clusters)
        cset$sorting <- sorting
        return(cset)
    }
    
    ## each clustering in \code{cset} comes with a cluster-cluster
    ## similarity matrix (used for scoring function \code{ccor});
    ## here we use it to get a rough sorting, simply starting with
    ## cluster '1'; the first cluster with the highest similarity
    ## to cluster '1' is taken as the next cluster
    for ( k in seq_len(ncol(cset$clusters)) ) {
        Ccc <- cset$Ccc[[k]]
        #if ( is.null(colnames(Ccc)) )
        #  colnames(Ccc) <- rownames(Ccc) <- seq_len(ncol(Ccc)
        sorting <- colnames(Ccc)
        remaining <- sorting
        cl <- remaining[1]
        remaining <- remaining[remaining!=cl]
        cln.srt <- cl
        ## start at first cluster
        while ( length(remaining) > 1 ) {
            clcor <- Ccc[cl,remaining,drop=FALSE]
            ## get first cluster with highest correlation
            new <- colnames(clcor)[which.max(clcor)] 
            if ( verb>0 ) cat(paste("\t", cl, ">", new,
                                  round(max(clcor),2), "\n"))
            cl <- new
            remaining <- remaining[remaining!=cl]
            cln.srt <- c(cln.srt, new)
        }
        ## add last
        cln.srt <- c(cln.srt, remaining)
        cset$sorting[[k]] <- cln.srt
    }
    names(cset$sorting) <- colnames(cset$clusters)
    
    cset
}

#' Parameters for \code{\link{segmentCluster.batch}}.
#' 
#' Generates the parameter list (\code{varySettings}) for
#' \code{\link{segmentCluster.batch}}, using defaults
#' for all parameters not passed.
#' @inheritParams segmentClusters
#' @return Returns a parameter settings structure that can be used
#' in the batch function \code{\link{segmentCluster.batch}}.
#'@export
setVarySettings <- function(E=c(1,3),
                            S="ccor",
                            M=100,
                            Mn=100,
                            a=-2, nui=c(1,3),
                            nextmax=TRUE,
                            multi="max",
                            multib="max") {
    list(E=E, S=S, M=M, Mn=Mn, a=a, nui=nui, # scoring
         nextmax=nextmax, multi=multi, multib=multib) # backtracing
}

#' Batch wrapper for \code{\link{segmentClusters}}.
#' 
#' A high-level wrapper for multiple runs of segmentation by
#' \code{\link{segmentClusters}} for multiple clusterings and/or
#' multiple segmentation parameters. It additionally allows to
#' tag adjacent segments to be potentially fused due to similarity
#' of their clusters.
#' 
#' @details This is a high-level wrapper for \code{\link{segmentClusters}}
#' which allows segmentation over multiple clusterings as provided by the
#' function \code{\link{clusterTimeseries}} and over multiple segmentation
#' parameters. Each parameter in the list \code{varySettings} can be
#' a vector and ALL combinations of the passed parameter values will
#' be used for one run of \code{\link{segmentClusters}}.
#' The resulting segment table, list item "segments" of the returned object,
#' is a \code{\link[base:data.frame]{data.frame}} with additional
#' columns "ID" and "type", automatically generated strings indicating
#' the used parameters (each "type" reflects one parameter set), and
#' "colors", indicating the automatically generated color of the assigned
#' cluster label.
#' @param cset a clustering set as returned by \code{\link{clusterTimeseries}}
#' @param varySettings list of settings where each entry can be a vector;
#' the function will construct a matrix of all possible combinations of
#' parameter values in this list, call \code{\link{segmentClusters}} for
#' each, and report a matrix of segments where the segment `type' is
#' constructed from the varied parameters; see option \code{short.name}.
#' A varySettings list with all required (default) parameters can be
#' obtained via function \code{\link{setVarySettings}}.
#' @param fuse.threshold if adjacent segments are associated with clusters
#' the centers of which have a Pearson correlation \code{>fuse.threshold}
#' the field "fuse" will be set to 1 for the second segments (top-to-bottom
#' as reported)
#' @param rm.nui remove nuisance cluster segments from final results
#' @param type.name vector of strings selecting the parameters which will be
#' used as segment types. Note, that all parameters that are actually varied
#' will be automatically added (if missing). The list can include parameters
#' from time-series processing found in the "clustering" object \code{cset}
#' as \code{cset$tids}.
#' @param short.name default type name construction; if TRUE (default)
#' parameters that are not varied will not be part of the segment type and ID.
#' This argument has no effect if argument \code{type.name} is set.
#' @param id if set, the default segment IDs, constructed from numbered
#' segment types, are replaced by this
#' @param save.matrix store the total score matrix \code{S(i,c)} and the
#' backtracing matrix \code{K(i,c)}; useful in testing stage or for
#' debugging or illustration of the algorithm;
#' TODO: save.matrix is currently not implemented, since batch function
#' returns a matrix only
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @return Returns an object of class "segments", just as its base function
#' \code{\link{segmentClusters}}, but the main segment table, list item
#' "segments", is a \code{\link[base:data.frame]{data.frame}} with additional
#' columns "ID" and "type", automatically generated strings indicating
#' the used parameters (each "type" reflects one parameter set), and
#' "colors", indicating the automatically generated color of the assigned
#' cluster label.
#' @examples
#' # load example data, an RNA-seq time-series data from a short genomic
#' # region of budding yeast
#' data(primseg436)
#' 
#' # 1) Fourier-transform time series:
#' tset <- processTimeseries(ts=tsd, na2zero=TRUE, use.fft=TRUE,
#'                           dft.range=1:7, dc.trafo="ash", use.snr=TRUE)
#'
#' # 2) cluster time-series several times into K=12 clusters:
#' cset <- clusterTimeseries(tset, K=c(12,12,12))
#'
#' # 3) choose parameter ranges, here only E is varied 
#' vary <- setVarySettings(M=100, E=c(1,3), nui=3, S="icor", Mn=20)
#' 
#' # 4) ... segment ALL using the batch function:
#' \dontrun{ ## NOTE: takes too long for CRAN example timing restrictions
#' segments <- segmentCluster.batch(cset=cset, varySettings=vary)
#' 
#' # 5) inspect results:
#' print(segments)
#' plotSegmentation(tset, cset, segments)
#' 
#' # 6) and get segment border table. Note that the table has
#' #    additional columns "ID" and "type", indicating the used parameters,
#' #    and "color" providing the color of the cluster the segment was
#' #    assigned to. This allows to track segments in the inspection plots.
#' sgtable <- segments$segments
#' }
#' 
#'@export
segmentCluster.batch <- function(cset, varySettings=setVarySettings(),
                                 fuse.threshold=0.2, rm.nui=TRUE, 
                                 type.name, short.name=TRUE, id,
                                 save.matrix=FALSE, verb=1) {

    ## TODO: allow defaults; getSettings to get full list!

    ## SETTING UP PARAMETER MATRIX
    ## 1) combine clusterings with segmentation parameters
    nk <- length(cset$K)
    vS <- append(list(K=colnames(cset$clusters)), varySettings)
    vL <- sapply(vS,length)
    names(vL) <- names(vS)
    rL <- c(1,vL)
    params <- as.data.frame(matrix(NA,ncol=length(vS),
                                 nrow=prod(sapply(vS,length))))
    colnames(params) <- names(vS)
    ## fill parameter matrix
    for ( j in seq_len(ncol(params)) ) 
        params[,j] <- rep(rep(vS[[j]],prod(rL[seq_len(j)])),
                          each=prod(rL[(j+2):length(rL)]))

    ## TODO: add time-series processing info
    if ( !is.null(cset$tsid) ) {
        cllst <- strsplit(cset$tsid, "_")
        ## get class ids
        clid <- unique(unlist(lapply(cllst, function(x) sub(":.*","",x))))
        ## fill data.frame
        cltab <- data.frame(matrix(NA, nrow=length(cset$tsid),
                                   ncol=length(clid)))
        colnames(cltab) <- clid
        for ( i in seq_along(cllst) ) {
            tmp <- strsplit(cllst[[i]],":")
            class <- unlist(lapply(tmp, function(x) x[2]))
            names(class) <- unlist(lapply(tmp, function(x) x[1]))
            cltab[i, names(class)] <- class
        }
        ## copy each to all existing param values
        if ( any(rL>1) ){ # do only if necessary, because otherwise
                          # cltab will be converted to a vector!
            cltab <- apply(cltab, 2, function(x)
                           rep(rep(x,prod(rL[1])),
                               each=prod(rL[(1+2):length(rL)])))
        }
        params <- cbind(cltab, params)
        ## recalculate vL; number of types in each setting
        vL <- apply(params,2,function(x) length(unique(x)))
        names(vL) <- colnames(params)
    }

    ## TODO: store previous processing IDs explicitly, to be used in plots

    ## CONSTRUCT SEGMENT CLASS NAMES!
    if ( missing(type.name) ) { 
        type.name <- colnames(params)
        ## rm those with length==1 to keep short names
        ## UNLESS there is no variation
        if ( short.name )
          if ( sum(vL>1)>0 )
            type.name <- type.name[vL>1]
          else
            type.name <- "S" # DEFAULT ID: scoring function
    } else {
        ## add all varied parameters appear in type.name
        varied <- names(vL[vL>1])
        type.name <- unique(c(type.name, varied))
    }
    
    if ( verb>0 )
        cat(paste("SEGMENTATIONS\t",nrow(params),"\n",sep=""))

    allsegs <- sgtypes <- NULL
    ## internal matrices (S, K, S1)
    if ( save.matrix ) 
      SK <- rep(list(NA), nrow(params))
    else
      SK <- NULL
    ## colors & sorting
    seg.col <- rep(list(NA), nrow(params))
    seg.srt <- rep(list(NA), nrow(params))
    ## run time
    elapsed <- rep(NA, nrow(params))
    
    ## TODO: convert this loop to lapply and try parallel use!
    ## TODO: redirect messages to msgfile or store in results
    for ( i in seq_len(nrow(params)) ) {

        ## construct segment type name  
        sgtype <- paste(paste(type.name,
                              unlist(lapply(params[i,type.name],as.character)),
                              sep=":"),
                        collapse="_")
        ## clustering K comes formatted already (K:id)
        sgtype <- sub("K:K:","K:", sgtype)
        sgtypes <- c(sgtypes, sgtype)

        ## clustering input
        K <- as.character(params[i,"K"])
        seq <- cset$clusters[,K]


        ## scoring params
        S <- as.character(params[i,"S"])
        E <-  as.numeric(params[i,"E"])
        M <-  as.numeric(params[i,"M"])
        Mn <- as.numeric(params[i,"Mn"])
        nui<- as.numeric(params[i,"nui"])
        a <-  as.numeric(params[i,"a"])

        ## back-tracing params
        multi   <- as.character(params[i,"multi"])
        multib  <- as.character(params[i,"multib"])
        nextmax <- as.logical(params[i,"nextmax"])

        if ( S=="ccor" ) csim <- cset$Ccc[[K]]
        if ( S=="icor" ) csim <- cset$Pci[[K]]
        if ( S=="ccls" ) csim <- NULL

        if ( verb>0 )
            cat(paste("SEGMENT TYPE\t",sgtype,
                      "\t", i,"of",nrow(params),"\n",sep=""))

        ## TODO: pass cset, inherit colors and set type ID there!
        seg <-segmentClusters(seq=seq,csim=csim,E=E,
                              S=S,M=M,Mn=Mn,nui=nui,a=a,
                              multi=multi,multib=multib,nextmax=nextmax,
                              save.matrix=save.matrix,rm.nui=rm.nui,verb=verb)

        ## retrieve algo-internal vectors and matrices
        ## S1: s(1,i,C) - scoring function from j=1 to i
        ## S: S(i,C) - score matrix
        ## K: K(i,C) - backtracing
        ## all can be plotted via dedicated functions
        if ( save.matrix ) 
            SK[[i]] <- seg$SK[[1]]

        ## pass on cluster coloring as segment coloring
        if ( "colors" %in% names(cset) )
          seg.col[[i]] <- cset$colors[[K]]
        if ( "sorting" %in% names(cset) )
          seg.srt[[i]] <- cset$sorting[[K]]

        ## run.time
        if ( "elapsed"%in%names(seg))
            elapsed[i] <- seg$elapsed
        
        ## tag adjacent segments from correlating clusters
        close <- fuseTagSegments(seg$segments, Ccc=cset$Ccc[[K]],
                                 fuse.threshold=fuse.threshold)
        if ( sum(close)>0 & verb>0 )
            cat(paste("Fused tags\t",sum(close), "\n",sep=""))

        ## collect results
        if ( nrow(seg$segments) > 0 ) {

            ## add colors as column
            ## TODO: redundant with cluster colors; but both are used
            colors <- seg$segments[,"CL"]
            if ( "colors" %in% names(cset) )
                colors <- cset$colors[[K]][as.character(seg$segments[,"CL"])]

            ##close <- rep(FALSE, nrow(seg$segments))
            sgids <- paste(sgtype, seq_len(nrow(seg$segments)),sep="_")
            segs <- data.frame(ID=sgids,
                               type=rep(sgtype,length(sgids)),
                               seg$segments,
                               fuse=close,
                               color=colors, stringsAsFactors=FALSE)
            allsegs <- rbind(allsegs,segs)
            
        } 
        if ( verb>0 )
            cat(paste("Segments\t", nrow(seg$segments), "\n",sep=""))
    }
    if ( is.null(allsegs) & verb>0 )
        cat(paste("Total segments\t0\n",sep=""))
    else {
        cat(paste("Total segments\t",nrow(allsegs),"\n",sep=""))
        ## OVERRIDE ID
        if ( !missing(id) ) 
            allsegs[,"ID"] <- paste(id, seq_len(nrow(allsegs)), sep="_")
    }

    ## name all stored data!
    if ( save.matrix )
      names(SK) <- sgtypes
    rownames(params) <- names(seg.col) <- names(seg.srt) <- sgtypes

    ## store length of sequence
    N <- nrow(cset$clusters)
    
    ## TODO: introduce and use classes for segment results
    sset <- list(segments=allsegs, N=N, colors=seg.col, sorting=seg.srt,
                 SK=SK, settings=params, ids=sgtypes, elapsed=elapsed)
    class(sset) <- "segments"
    return(sset)
}

## tags adjacent segments if they are from correlating (>\code{fuse.thresh})
## clusters
fuseTagSegments <- function(segs, Ccc, fuse.threshold=.2) {

    if ( nrow(segs)==0 ) return(NULL)
    if ( nrow(segs)==1 ) return(FALSE)
    
    fuse <- rep(-Inf,nrow(segs))
    for ( j in 2:nrow(segs) ) {
        previous <- segs[j-1,1]
        current <- segs[j,1]
        if ( !0 %in% c(previous,current) ) 
            fuse[j] <- Ccc[current,previous]
    }
            
    ## FUSE directly adjacent if clusters correlate?
    adj <- segs[2:nrow(segs),2] - segs[2:nrow(segs)-1,3] ==1
    close <- c(FALSE,adj) & fuse > fuse.threshold

    close
}

