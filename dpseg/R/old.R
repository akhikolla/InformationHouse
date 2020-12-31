### OLD FUNCTION OF FIRST IMPLEMENTATIONS,
### NOT USED BUT IMPLEMENTATION WAS MORE GENERAL & MODULAR
### PERHAPS FOR LATER RE-USE & EDUCATIONAL PURPOSES

## linear regression function, RcppEigen::fastLm,
## must return an S3 object containing the residuals,
## intercept and slope as coef(linregf())[1:2] are currently
## only required to add slopes, and in the current predict function
linregf <- function(x, y) 
    lm(y ~ x) 

## fitscore function for linear regression, 
## fitscore: -variance of residuals of linear fit from i+1 to j
varscore <- function(x, y, fitf=linregf, ...) 
    -var(fitf(x, y, ...)$residuals)

## calculate fitscore matrix
fitscore <- function(x, y, minl=0, maxl=Inf,
                     fitf=linregf, scoref=varscore) {
    
    ## TODO: switch to apply and use parallelization?
    ## ... or move to Rcpp alltogether

    EPS <- matrix(NA, nrow=length(x), ncol=length(x))
    for ( j in seq_along(x) )
        ## TODO: fix this inner for loop!
        if ( j-minl >=  max(1,j-maxl) )
            for ( i in max(1,j-maxl):(j-minl) ){
                rng <- (i+1):j
                yn <- y[rng]
                xn <- x[rng]
                ## remove NA or infinite data (used on log of exp. data!)
                xn <- xn[!is.na(yn)&is.finite(yn)]
                yn <- yn[!is.na(yn)&is.finite(yn)]
                if ( length(yn)>1 ) 
                    EPS[i,j] <- scoref(xn, yn, fitf=fitf)
            }
    invisible(EPS)

}

## direct calculation of error function after PFS, 20181219
## this was then used as a basis for the new Rcpp implementation
fitscored <- function(x, y, minl=0, maxl=Inf, scoref, fitf) {

    vect <- rep(NA, length(x))
    sy <- sx <- sy2 <- sx2 <- sxy <- vect
    sy2[1] <- y[1]^2
    sx2[1] <- x[1]^2
    sy[1] <- y[1]
    sx[1] <- x[1]
    sxy[1] <- y[1]*x[1]
    eps <- matrix(NA, nrow=length(x), ncol=length(x)) # former FS!
    for ( j in 2:length(x) ) {
        sy2[j] <- sy2[j-1] + y[j]^2      # \sum y^2
        sy[j] <-  sy[j-1] + y[j] 
        sx2[j] <- sx2[j-1] + x[j]^2      # \sum x^2
        sx[j] <-  sx[j-1] + x[j] 
        sxy[j] <- sxy[j-1] + x[j]*y[j]   # \sum x*y
        if ( j-minl >=  max(1,j-maxl) )
            for ( i in max(1,j-maxl):(j-minl) ){
                n <- j-i
                ## Syy = \sum y^2 - (\sum y)^2/n
                Syy <- (sy2[j] - sy2[i]) - (sy[j]-sy[i])^2/n
                ## Sxx = \sum x^2 - (\sum x)^2/n
                Sxx <- (sx2[j] - sx2[i]) - (sx[j]-sx[i])^2/n
                ## Sxy = \sum x*y - sum(x)*sum(y)/n
                Sxy <- (sxy[j] - sxy[i]) - (sx[j]-sx[i])*(sy[j]-sy[i])/n
                ## variance, from i to j:
                eps[i,j] <- (Syy - Sxy^2/Sxx) / (n-1)
            }
    }
    invisible(-eps)
}


## fitscore version using apply, BUT NOT REALLY FASTER!
## was tested in dpseg_example.R
fitscorel <- function(x, y, minl=0, maxl=Inf,
                      fitf=linregf, scoref=varscore) {
    ## calculate fits using apply loops
    fsl <- lapply(seq_along(x), function(j) {
        fs <- rep(NA, length(x))
        if ( j-minl >=  max(1,j-maxl) ) {
            rng <-max(1,j-maxl):(j-minl) 
            fs[rng] <- sapply(rng, function(i) {
                rng <- (i+1):j
                yn <- y[rng]
                xn <- x[rng]
                ## remove NA or infinite data (used on log of exp. data!)
                xn <- xn[!is.na(yn)&is.finite(yn)]
                yn <- yn[!is.na(yn)&is.finite(yn)]
                if ( length(yn)>1 ) 
                    return(scoref(xn, yn, fitf=fitf))
                else (return(NA))
            })
        }
        fs
    })
    EPS <- do.call(cbind, fsl)
    invisible(EPS)
}

## recursion: calculate S_j
## 
## The recursion calculates
## \code{S_j = max(i<(j-minl) S_i + fitscore(i+1,j) - P},
## where \code{P} is a jump penality that implicitly regulates the
## number of segments, \code{minl} and \code{maxl} are minimal and
## maximal lengths of segments.
## @export
score <- function(EPS, minl=3, maxl, P=0) {

    ## RECURSION S_j = max(i<(j-minl) S_i + fitscore(i+1,j) - P

    N <- ncol(EPS)

    ## maximal length?
    if ( missing(maxl) ) maxl <- N
    
    ## initialize
    S <- imax <- rep(0, N)
    S[seq_len(minl)] <- 0 - P

    ## calculate score S_j
    for ( j in (minl+1):N ) {

        ## calculate all scores i->j
        si <- rep(-Inf, j-minl)
        for ( i in max(1,j-maxl):(j-minl) ) 
            si[i] <- S[i] + EPS[i,j] - P
        
        ## get maximum 
        S[j] <- max(si)
        ## ... and store which i yielded maximum
        imax[j] <- which.max(si)
        
    }
    
    invisible(list(imax=imax, S=S))
}

traceback_old <- function(S) {

    imax <- S[["imax"]]
    
    ## backtracing
    end <- imax[length(imax)]
    ends <- c(end,length(S$S))
    while( end>0 ) {
        end <- imax[end]
        ends <- c(end, ends)
    }

    ## convert to start:end coordinates
    ## imax[j] for max(i<(j-minl) Si + fitscore(i+1,j)
    segs <- cbind(start=ends[2:length(ends)-1]+1,
                  end=ends[2:length(ends)])
    invisible(segs)
}

## get lm fits for final segments
## to be used also for generic predict
fitseg_old <- function(x,y,segments, fitf=linregf) {

    ## TODO: use try to catch and store official errors!
    fits <- rep(list("fiterror"), nrow(segments))

    for ( i in seq_len(nrow(segments)) ) {
        rng <- segments[i,"start"]:segments[i,"end"]
        yn <- y[rng]
        xn <- x[rng]
        ## remove NA or infinite data (used on log of exp. data!)
        xn <- xn[!is.na(yn)&is.finite(yn)]
        yn <- yn[!is.na(yn)&is.finite(yn)]
        if ( length(yn)>1 ) 
            fits[[i]] <- fitf(xn, yn) # defined herein!
    }
    fits
}

#' inefficient \code{\link{dpseg}} implementation
#'
#' See  \code{\link{dpseg}} for a current version of this algorithm.
#' Note: this was a first test implementation of the linear
#' piecewise segmentation by a dynamic programming approach.
#' This implementation is very slow. A much more efficient
#' version, \code{\link{dpseg}}, calculates the variance of residuals of a
#' linear regression incrementally while looping through the recursion, and
#' is implemented in \code{Rcpp}.
#' See there for details on the algorithm.
#' This version is kept alive, since it is a more general implementation,
#' allowing to test different regression and scoring functions by
#' command-line arguments.
#'
#' The recursion calculates \eqn{S_j = max ( S_i +
#' fitscore(i+1,j)) - P}, where the fitscore is the variance of the
#' residuals of a linear regression (\code{lm(y~x})) between
#' \eqn{x_{i+1}} to \eqn{x_j}, \code{P} is a jump penality that
#' implicitly regulates the number of segments, \code{minl} and
#' \code{maxl} are minimal and maximal lengths of segments. Uses
#' \code{\link[RcppEigen:fastLm]{RcppEigen:fastLm}} for linear
#' regression.
#' @param x x-values
#' @param y y-values
#' @param minl minimal segment length
#' @param maxl maximal segment length
#' @param P jump penalty, increase to get fewer segments #
#'     @inheritParams score
#' @param EPS a pre-calculated fitscore matrix, will be generated if
#'     missing
#' @param store.matrix store the fitscore matrix
#' @param fitscoref the heavy-load loop that fills the fitscore matrix
#'     using \code{fitf} and \code{scoref}
#' @param fitf fit function, used in the scoring function
#'     \code{scoref}; (TODO: currently expecting a fit object that
#'     provides intercept and slope as \code{coef(obj)[1:2]} only for
#'     the result table)
#' @param scoref function to calculate a score from the passed fit
#'     function
#' @param verb print progress messages
#' @return Returns a list of result structures very similar to the
#'     list of class "dpseg" returned by function \code{\link{dpseg}},
#'     except for the name of the scoring function matrix, here:
#'     \code{EPS}. See \code{?dpseg} for detailed information on these
#'     structures.
#' @examples
#' \donttest{
#'
#' ## NOTE: not run because it's too slow for R CMD check --as-cran
#' ## calculate linear segments in semi-log bacterial growth data
#' ## NOTE: library loads bacterial growth curve data as data.frame oddata
#' Sj <- dpseg_old(x=oddata$Time, y=log(oddata$A3), minl=5, P=0.0001, verb=1)
#' 
#' ## inspect resulting segments
#' print(Sj)
#' 
#' ## plot results
#' plot(Sj, delog=TRUE, log="y")
#' 
#' ## NOTE: predict method & movie function do not work for dpseg_old
#' }
#' @export
dpseg_old <- function(x, y, minl, maxl=length(x), P=0, EPS, store.matrix=FALSE,
                     fitscoref=fitscore, fitf=linregf, scoref=varscore,
                     verb=0) {

    ## 1) construct fitscore matrix
    if ( missing(EPS) ) {
        if ( verb>0 )
            cat(paste("pre-calculating fit score matrix; can take a while\n"))
        EPS <- fitscoref(x,y, minl=minl, maxl=maxl,
                        scoref=scoref, fitf=fitf)
    }
    ## 2) recursion
    ## RECURSION: S_j = max(i<(j-minl) S_i + fitscore(i+1,j) - JUMPPENALTY
    if ( verb>0 )
        cat(paste("calculating recursion\n"))
    Sj <- score(EPS, minl=minl, maxl=maxl, P=P)

    ## 3) back-tracing
    if ( verb>0 )
        cat(paste("back-tracing\n"))
    segs <- traceback_old(Sj)

    ## 4) add fits for final segments
    fits <- fitseg_old(x, y, segs, fitf=fitf)
    
    ## 5) add slopes - GROWTHRATES, but
    ## TODO: make this more general, since slopes
    ## are only available for linear regressions
    segs <- slopes(x, y, segs, fits)

    ## add matrix and segments to results
    Sj[["EPS"]] <- EPS
    Sj[["segments"]] <- segs

    ## add original data, parameters and segment fits
    Sj[["xy"]] <- grDevices::xy.coords(x=x,y=y)
    Sj[["parameters"]] <- c(minl=minl,maxl=maxl,P=P)
    Sj[["fits"]] <- fits
    
    ## re-order
    ord <- c("segments","xy","parameters","S","imax","fits")
    if ( store.matrix ) 
        ord <- c("segments","xy","parameters","S","imax","EPS", "fits")
        
    Sj <- Sj[ord]

    class(Sj) <- "dpseg"
    invisible(Sj)
}


## UNUSED OLD R IMPLEMENTATION OF THE DYN.PROG. RECURSION
## NOW IN src/segment.cpp; can be used by passing dpseg:::scoref
## as argument scorefnc to dpseg
## NOTE: this ignores minl at the beginning!
scoref_oldr <- function(x, y, minl, maxl, P, storem=FALSE) {
    
    N <- length(x)

    ## A) INITIALIZATION
    ## initialize Scoring vector
    S <- imax <- ieps <- rep(0, N)
    S[seq_len(minl)] <- 0 - P
    
    ## initialize x2, y2, xy
    sy <- sx <- sy2 <- sx2 <- sxy <- rep(NA, N)
    sx[1] <- x[1]
    sy[1] <- y[1]
    sx2[1] <- x[1]^2
    sy2[1] <- y[1]^2
    sxy[1] <- y[1]*x[1]
    
    ## store matrix: only required for educational movies!
    if ( storem ) 
        EPS <- matrix(NA, nrow=length(x), ncol=length(x)) # former FS!
    
    
    ## B) RECURSION
    for ( j in 2:length(x) ) {
        
        ## fill up variances
        sx[j] <-  sx[j-1]  + x[j]        # \sum x
        sy[j] <-  sy[j-1]  + y[j]        # \sum y
        sx2[j] <- sx2[j-1] + x[j]^2      # \sum x^2
        sy2[j] <- sy2[j-1] + y[j]^2      # \sum y^2
        sxy[j] <- sxy[j-1] + x[j]*y[j]   # \sum x*y
        
        ## only calculate S[i] if required
        if ( j-minl >=  max(1,j-maxl) ) {
            ## calculate all scores i->j
            si <- rep(-Inf, j-minl)
            epsi <- rep(0, j-minl)

            ## analyze Sxx==0 for yeast chr2, j=476359, i=476354
            ##for (j in 476359 )

            for ( i in max(1,j-maxl):(j-minl) ){
                n <- j-i
                ## Syy = \sum y^2 - (\sum y)^2/n
                Syy <- (sy2[j]-sy2[i]) - (sy[j]-sy[i])^2/n
                ## Sxx = \sum x^2 - (\sum x)^2/n
                Sxx <- (sx2[j]-sx2[i]) - (sx[j]-sx[i])^2/n
                ## Sxy = \sum x*y - sum(x)*sum(y)/n
                Sxy <- (sxy[j]-sxy[i]) - (sx[j]-sx[i])*(sy[j]-sy[i])/n

                ##cat(paste(Sxx,"\n"))
                ##if ( Sxx==0 ) stop(i)
                ##}
                
                ## SCORING FUNCTION: variance, from i to j:
                epsi[i] <- (Syy - Sxy^2/Sxx) / (n-1)
                
                ## calculate score
                if ( Sxx!=0 ) # TODO: why can Sxx be 0?
                    si[i] <- S[i] - epsi[i] - P
                
                if ( storem )
                    EPS[i,j] <- -epsi[i]
            }
            
            ## get maximum 
            S[j] <- max(si)
            ## ... and store which i yielded maximum
            imax[j] <- which.max(si)
            ## ... and variances i+1 to j
            ieps[j] <- epsi[imax[j]]
        }
     }

    if ( storem ) invisible(list(S=S, imax=imax, ivar=ieps, EPS=EPS))
    else invisible(list(S=S, imax=imax, values=list(ivar=ieps)))
}
