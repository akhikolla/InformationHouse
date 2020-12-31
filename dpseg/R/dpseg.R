#' dpseg : linear segmentation by dynamic programming
#' @docType package
#' @name dpseg
#' @author Rainer Machne \email{machne@hhu.de}, Peter F. Stadler \email{studla@bioinf.uni-leipzig.de}
#' 
## @description Solves
## \code{S_j = max_{j-maxl < i <= j-minl} S_i + score(i+1,j) - P}, with
## score: -variance of residuals of a linear regression from i+1 to j, and
## P: a jump penality that implicitly regulates length.
#' @section Dependencies: The package strictly depends only on \code{RcppEigen}.
#' All other dependencies are usually present in a
#' basic installation (\code{stats}, \code{graphics}, \code{grDevices}).
## @importFrom RcppEigen fastLm
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline arrows axis lines mtext par plot points
#' @importFrom stats lm coef var predict cor residuals smooth.spline median
#' @importFrom Rcpp evalCpp 
#'@useDynLib dpseg
NULL # this just ends the global package documentation



#' backtracing \code{dpseg} segment break-points
#'
#' Backtracing segment borders from the \code{imax} vector of a
#' \code{\link{dpseg}} recursion. This function is implemented more
#' efficiently in Rcpp; the R code is kept for documentation,
#' benchmarking and development.
#' @param imax integer vector of segment borders as returned by
#'     \code{dpseg} recursion functions
#' @param jumps allwo discontinuous jumps: move 1 index position back,
#'     only for \eqn{S_{i-1} + score(i,j)}
#' @return an integer vector of segment ends
#' @export
backtrace_r <- function(imax, jumps=0) {
    
    end <-length(imax) # start at j==N
    ends <- c(end) 
    while( end>1 ) {
        ## get i that was max S[i - jumps] score(,j)
        end <- imax[end] - jumps
        ends <- c(end, ends)
    }
    ends
}

#' construct a segment table
#'
#' Constructs a segment table from segment ends (\code{imax}) returned
#' by \code{\link{dpseg}} backtracing functions
#' \code{\link{backtrace_r}} and \code{backtrace_c}. Correct segment
#' break-points require to know whether segment recursion was run with
#' the \code{jumps} option of \code{\link{dpseg}}. In joint segments
#' (\code{jumps=FALSE}) segment borders are part of both left and
#' right segments.
#' @param ends integer vector of segment ends
#' @param starts integer vector of segment starts
#' @param jumps same parameter as passed to recursion function,
#'     allowing for discontinous jumps (TRUE) or enforcing joint
#'     segments (FALSE)
#' @return a table with segment \code{start} and \code{end} columns
#' @export
sgtable <- function(ends, starts, jumps=TRUE) {

    ## ends/starts are the vectors returned by backtrace
    ## note that this vector consists of imax breakpoints + first/last index
    if ( missing(starts) ) {
        if ( length(ends)==2 ) { # short segment
            starts <- ends[1]
            ends <- ends[2]
        } else {
            if ( jumps ) # S[i-1] score(i,j) or S[i] score(i+1,j)
                starts <- c(ends[1], ends[3:length(ends)-1]+1)
            else  # S[i] score(i,j)
                starts <- c(ends[1], ends[3:length(ends)-1])
            ends <- ends[2:length(ends)]
        }
    }
    if ( missing(ends) ) {
        if ( length(starts)==2 ) {
            ends <- starts[2]
            starts <- starts[1]
         } else {
            if ( jumps )
                ends <- c(starts[3:length(starts)-1]-1,
                          starts[length(starts)])
            else
                ends <-  c(starts[3:length(starts)-1],
                           starts[length(starts)])
            starts <- starts[2:length(starts)-1]
        }
    }
    cbind(start=starts, end=ends)
}


#' dpseg : linear segmentation by dynamic programming
#' 
#' \code{dpseg} splits a curve (x,y data) into linear segments by a
#' straight forward dynamic programming recursion: \deqn{S_j =
#' max(S_{i-jumps} + score(i,j) -P)} where score is a measure of the
#' goodnes of the fit of a linear regression (equiv. to
#' \code{lm(y~x)}) between data points \eqn{i<j}. The default scoring
#' function is simply the negative variance of residuals of the linear
#' regression (see arguments \code{type} and \code{scoref}). \code{P}
#' is a break-point penality that implicitly regulates the number of
#' segments (higher P: longer segments), and \code{jumps==1} allows
#' for disjoint segments.  The arguments \code{minl} and \code{maxl}
#' specify minimal (\eqn{i\le j-minl}) and maximal (\eqn{i\ge j-maxl})
#' segment lengths, which allows to significantly decrease memory
#' usage when expected segment lengths are known.
#'
#' See the \code{vignette("dpseg")} for the theory and details on the
#' choice of scoring functions and selection of the penalty parameter
#' P.
#' @param x x-values, not used if \code{y} is a scoring function
#'     matrix
#' @param y y-values, or a pre-calculated scoring function matrix
#'     \eqn{SCR_{i,j}} (eg. from a previous run of \code{dpseg}). See
#'     section "Value" below for details on the structure \eqn{SCR_{i,j}}.
#' @param maxl maximal segment length, \eqn{i\ge j-maxl}
#' @param minl minimal segment length, \eqn{i\le j-minl}
#' @param jumps allow for jumps between segments, if \code{TRUE}
#'     segment ends are 1 index left of the segment starts
#' @param P break-point penalty, increase to get longer segments with
#'     lower scores (eg. higher residual variance)
#' @param S0 initialization of \eqn{S_0}, choose high enough to avoid
#'     length 1 cutoffs at start
#' @param type type of scoring function: available are "var" for
#'     "variance of residuals", "cor" for Pearson correlation, or "r2"
#'     for r-squared; see the package \code{vignette("dpseg")} for
#'     details.
#' @param scoref alternative scoring function
#' @param store.values store scoring values (linear regression
#'     results)
#' @param store.matrix store the fitscore matrix
#' @param verb print progress messages
#' @param add.lm add a linear fit using R base \code{lm} for final
#'     segments; may save memory/speed if \code{store.values==FALSE}
#' @param move logical indicating whether move is required in
#'     backtracing, required for the alternative recursion \eqn{S_i +
#'     score(i+1,j)}
#' @param recursion internal recursion function to be used for
#'     segmentation; used for debugging, benchmarking and development,
#'     and required for putative novel scoring functions \code{scoref}
#' @param backtrace internal function to be used for back-tracing;
#'     used for debugging, benchmarking and development, and may be
#'     required to test novel scoring functions \code{scoref} and/or
#'     \code{recursion}
#' @param ... further arguments to \code{recursion}
#' @return Returns a list object of class \code{dpseg} (with
#'     \code{print.dpseg} \code{plot.dpseg} and \code{predict.dpseg}
#'     methods). The main result of the algorithm is a table
#'     (\code{data.frame}) of predicted segments in list object
#'     \code{segments}. The original data, run parameters and
#'     (optionally) additional data calculated and used by the
#'     algorithm are also returned.
#'
#' \describe{
#'
#' \item{segments:}{main result table: a \code{data.frame} that lists
#'     the start and end x-values of the segments, the start and end
#'     indices (i,j) in the data vectors, the linear regression
#'     coefficients and goodness-of-fit measures for the segments
#'     (intercept, slope, r-squared, variance of residuals). If
#'     \code{dpseg} was called with a pre-calculated scoring matrix,
#'     the table only contains start and end indices i,j. If option
#'     \code{add.lm=TRUE} or the result object was sent through
#'     function \code{\link{addLm}} the table additionally contains
#'     results from R's \code{lm}, indicated by an ".lm" suffix.}
#'
#' \item{S:}{results of the recursion, ie. \eqn{S_j} in above
#'     equation.}
#'
#' \item{imax:}{vector \eqn{j=1,\dots,n}, storing the \eqn{i_{max}}
#'     that yielded \eqn{S_j}, ie., the sole input for the backtracing
#'     function.}
#'
#' \item{values:}{linear regression coefficients and measures for the
#'     segment ending at \eqn{j} and starting at
#'     \eqn{i_{max}(j)}. Only present if \code{store.valus=TRUE}.}
#'
#'
#' \item{SCR:}{scoring function matrix \eqn{SCR_{i,j} = score(i,j)}
#'     where positions j are the columns and i the rows; a banded
#'     matrix with non-NA values between \eqn{i\le j-minl} and
#'     \eqn{i\ge j-maxl}. Note, that this matrix can be re-used in
#'     subsequent calls as \code{dpseg(y=previous$SCR)} which runs
#'     much faster and allows to efficiently scan for alternative
#'     parameters. Only present if \code{store.matrix=TRUE}.}
#'
#' \item{fits:}{result objects from \code{lm}. Only present if
#'     \code{add.lm=TRUE}.}
#'
#' \item{traceback:}{result of the call to the backtracing function:
#'     ends of the segments.}
#'
#' \item{xy:}{original x/y data (xy.coords).}
#'
#' \item{removed:}{index of NA/Inf values that were removed before
#'     running the alorithm.}
#'
#' \item{parameters:}{used parameters \code{P}, \code{jumps},
#'     \code{maxl} and \code{minl}.}
#' }
#' @examples
#' 
#' ## calculate linear segments in semi-log bacterial growth data
#' ## NOTE: library loads bacterial growth curve data as data.frame oddata
#' segs <- dpseg(x=oddata$Time, y=log(oddata$A3), minl=5, P=0.0001, verb=1)
#' 
#' ## inspect resulting segments
#' print(segs)
#' 
#' ## plot results (also see the movie method)
#' plot(segs, delog=TRUE, log="y")
#' 
#' ## predict method
#' plot(predict(segs), type="l")
#' 
#' @export
dpseg <- function(x, y, maxl, jumps=FALSE, P=0, minl=3, S0=1,
                  type="var", scoref, verb=1, move, 
                  store.values=TRUE, store.matrix=FALSE, add.lm=FALSE,
                  recursion, backtrace, ...) {


    ## SELECT RECURSION
    
    ## new scoring function: use generic recursion, which expects a function
    ## of exactly type: scoref(x,y,i,j)
    if ( !missing(scoref) ) {

        fname <- as.character(substitute(scoref))
        if ( verb>0 )
            cat(paste0("using generic recursion for scoring function \"",
                       fname,"\"\n"))

        if ( missing(recursion) )
            recursion <- recursion_generic
        ## type is only used for default recursions!
        type <- fname
    }

    ## pre-calculated matrix for generic recursion
    ## TODO: instead introduce proper class for scoring function matrix 
    if ( "matrix" %in% class(y) ) { 
        
        if ( verb>0 )
            cat(paste("setup use with pre-calculated scoring matrix\n"))

        if ( missing(x) ) x <- seq_len(nrow(y))
        x.orig <- y.orig <- x  # dummy x
        rm <- rep(FALSE, length(x))

        ## required functions and parameters
        store.values <- store.matrix <- add.lm <- FALSE
        if ( missing(recursion) )
            recursion <- recursion_matrix # generic recursion function
        type <- "matrix"
    }

    ## 1) INPUT DATA PRE-PROCESSING

    ## x-y data
    ## rm NA and infinite values (may appear due to frequent log(y) input
    if ( !"matrix" %in% class(y) ) {
        if ( missing(x) ) x <- seq_along(y) # TODO: efficient equispaced algo
        x.orig <- x
        y.orig <- y
        rm <- is.na(y)|is.infinite(y)
        if ( sum(rm)>0 ) {
            if (  verb>0 )
                cat(paste("removing", sum(rm), " NA or infinite y value(s)\n"))
            warning("removing ", sum(rm), " NA or infinite y value(s)")
        }
        x <- x[!rm]
        y <- y[!rm]
        if ( length(y)<2 )
            stop("no finite non-NA y values left!")
    }

    ## DEFAULT SCORING AND BACKTRACE FUNCTIONS
    if ( missing(backtrace) )
        backtrace <- backtrace_c
    if ( missing(recursion) )
        recursion <- recursion_linreg_c


    ## ALLOW DISCONTINUOUS JUMPS?
    ## # only S[i] + score(i+1,j) requires jumps & !move
    ## TODO: move vs. jump was for testing different indexing,
    ##       keep or remove?
    if ( missing(move) ) move <- jumps 
    jumps <- as.numeric(jumps) # switch between (i,j) and (i+,j)
    move <- as.numeric(move) # switch between S[i-1] and S[i]
    
    ## total length
    N <- length(x)

    ## start at x=0 to avoid big integer overflow (x^2)
    ## and catastrophic cancellation
    x1 <- x[1]
    x <- x-x1
    
    ## correct maximal length
    if ( missing(maxl) ) maxl <- N - jumps + move
    if ( maxl>(N-jumps+move) ) maxl <- N - jumps + move
  
    ## correct minimal length
    if ( minl<3 ) {
        warning("minimal length ", minl, " too short, set to minl=3; N=",N)
        minl <- 3
    }
    if ( minl>=maxl ) {
        warning("minimal length ", minl, " too long, set to maxl-1; maxl=",maxl)
        minl <- maxl - 1
    }

    if ( sum(!is.na(y)) <= minl ) {
        warning("not enough data, returning NULL")
        return(NULL)
    }
    
    ## 2) DYNAMIC PROGRAMMING RECURSION
    parms <- data.frame(type=type, minl=minl, maxl=maxl,
                        P=P, jumps=jumps, stringsAsFactors=FALSE)
    if ( verb>0 )
        cat(paste("calculating recursion for", N, "datapoints\n"))
    if ( verb>1 )
        print(parms)

    if ( !missing(scoref) )
        Sj <- recursion(x=x, y=y,maxl=maxl, jumps=jumps, P=P, minl=minl, S0=S0,
                        type=type, storem=store.matrix, storev=store.values,
                        scoref=scoref, ...)
    else 
        Sj <- recursion(x=x, y=y,maxl=maxl, jumps=jumps, P=P, minl=minl, S0=S0,
                        type=type,storem=store.matrix,storev=store.values, ...)
        
    
    ## 2) BACKTRACING
    if ( verb>1 )
        cat(paste("done with recursion, back-tracing\n"))

    ## move=1 shifts ends from S[i-1] score(i,j)
    ends <- backtrace(Sj$imax, jumps=move)

    ## fix first segment in S[i-1] score(i,j)
    if ( move ) ends[1] <- ends[1]+1
    
    ## 3) RE-CREATE ORIGINAL DATA (removed Inf/NA values)

    ## TODO: generalize objects returned by recursion and backtrace
    ## (assign classes "backtrace", "values") to allow generic
    ## handling of removed NA/Inf values and passing of to backtrace function

    ## add back first x
    x <- x+x1

    ## add back NA/infinite values
    ## NOTE: makes dpseg very slow for many NA/Inf in long series!!
    ## TODO: find more efficient version or move to Rcpp 
    addValue <- function(x, r, val, add=0) {
        if ( r>length(x) )
            x <- c(x+add, val)
        else {
            x[seq(r,length(x))+1]  <- x[seq(r,length(x))] +add
            x[r] <- val
        }
        x
    }
    if ( any(rm) ) {
        if ( verb>0 ) cat(paste("adding back removed data points\n"))
        ridx <- which(rm)
        for ( i in seq_along(ridx) ) {
            r <- ridx[i]

            ## add to main recursion result
            Sj$S <- addValue(Sj$S, r, -Inf, add=0)
            Sj$imax <- addValue(Sj$imax, r, NA, add=1) # note: add 1
            ## fix segment ends
            ends[ends>=r] <- ends[ends>=r] +1

            ## add to all values
            if ( "values" %in% names(Sj) ) # new
                for ( id in names(Sj$values) )
                    Sj$values[[id]] <- addValue(Sj$values[[id]], r, NA, add=0)
            ## add to scoring matrix
            if ( store.matrix ) {
                SCR <- matrix(NA, nrow=length(Sj$S), ncol=length(Sj$S))
                SCR[-ridx,-ridx] <- Sj$SCR
                Sj$SCR <- SCR
            }
        }
        x <- x.orig
        y <- y.orig
    }

    if ( verb>1 )
        cat(paste("done with back-tracing, generating table\n"))

    ## 4) SEGMENT TABLE
    ## convert to start:end coordinates
    ## imax[j] for max((j-maxl)<i<=(j-minl) Si + fitscore(i+1,j)
    segs <- sgtable(ends, jumps=jumps)

    ## store parameters and matrices calculated during recursion
    if ( store.values & "values"%in%names(Sj) ) {

        val <- Sj$values
        ## add parameters

        ## TODO: get correct coordinates for first/last segment
        cor <- "end"
        maxi <- segs[,cor] 
        segs <- cbind(segs, intercept=val$icpt[maxi])
        segs <- cbind(segs, slope=val$islp[maxi])
        segs <- cbind(segs, r2=val$irsq[maxi])
        segs <- cbind(segs, var=val$ivar[maxi])

        ## TODO: perhaps remove the alternative scoring function and skip this?
        ## for S[i] + score(i+1,j)
        if ( "jvar"%in%names(val) ) {
            segs[1,"var"] <- val$jvar[segs[1,"end"]]
            segs[1,"intercept"] <- val$jcpt[segs[1,"end"]]
            segs[1,"slope"] <- val$jslp[segs[1,"end"]]
            segs[1,"r2"] <- val$jrsq[segs[1,"end"]]
        }
        ## correct intercepts!
        segs[,"intercept"] <- segs[,"intercept"] - segs[,"slope"]*x1
    } else if ( store.values )
        warning("recursion didn't return fit values: can't `store.values`")
    
    ## 5) optionally add R base lm fits, for debugging,
    ##    and for custom scoring functions/recursion_generic
    if ( add.lm ) {
        if ( verb>1 )
            cat(paste("done with tables, adding lm results\n"))

        ## perform lm fits
        fits <- fitseg(x, y, segs)
        ## add slopes 
        segs <- slopes(x, y, segs, fits)

        ## add full fit objects
        Sj[["fits"]] <- fits
    }
    
    ## add breakpoints in x coordinates!
    segs <- cbind.data.frame(x1=x[segs[,"start"]],
                             x2=x[segs[,"end"]],segs)

    ## 6) COLLECT RESULTS

    ## add matrix and segments to results
    Sj[["traceback"]] <- ends
    Sj[["segments"]] <- segs # main result table

    ## add original data, parameters and segment fits
    if ( "matrix" %in% class(y) ) Sj[["xy"]] <- list(x=x, y=y)
    else  Sj[["xy"]] <- grDevices::xy.coords(x=x, y=y)
    Sj[["removed"]] <- which(rm)
    parms <- data.frame(type=type, minl=minl, maxl=maxl,
                        P=P, jumps=jumps, stringsAsFactors=FALSE)
    Sj[["parameters"]] <- parms
    
    
    class(Sj) <- "dpseg"
    Sj
}

#' Adds linear regression data to \code{dpseg}
#' results or a table of segment borders.
#'
#' \code{addLm} takes a segment table (with start/end columns) or a
#' result object from code \code{\link{dpseg}}, calls base R function
#' \code{\link{lm}} for each segment, and adds slope, intercept, r2 and
#' variance of residuals to the segment table. This data is required
#' for plot and predict method, eg.  when \code{dpseg} was called with
#' a pre-calculated scoring matrix, or alternative scoring functions
#' or recursion.
#' @param dpseg result object (class "dpseg") returned by function
#'     \code{\link{dpseg}} or simply a segment table with "start" and
#'     "end" indices
#' @param x original x-data used
#' @param y original y-data used
#' @return Returns the input \code{dpseg} object or segment table, but
#'     with original \code{xy} data and \code{fit} results from a
#'     linear regression with base R (\code{lm(y~x)}) added to the
#'     results and linear regression coefficient and goodness of fit
#'     meaurs in the main \code{segments} table.
#' @examples
#'
#' ## 1: run dpseg with store.matrix=TRUE to allow re-rung
#' segs <- dpseg(x=oddata$Time, y=log(oddata$A3), store.matrix=TRUE)
#'
#' ## 2: run dpseg with score function matrix input
#' segr <- dpseg(y=segs$SCR,  P=0.0001, verb=1)
#'
#' ## NOTE: only data indices i and j are provided in results
#' print(segr)
#' 
#' ## 3: add original data and linear regression for segments
#' ## NOTE: now also plot and predict methods work
#' segr <- addLm(segr, x=oddata$Time, y=log(oddata$A3))
#' print(segr)
#' 
#'@export
addLm <- function(dpseg,x,y) {
    
    if ( "dpseg" %in% class(dpseg) )
        segments <- dpseg$segments
    else segments  <- dpseg

    ## get or add xy!
    if ( missing(x) ) dpseg$xy$x
    if ( missing(y) ) dpseg$xy$y
    else {
        dpseg$xy <- grDevices::xy.coords(x=x, y=y)
        segments$x1 <- x[segments$start]
        segments$x2 <- x[segments$end]
    }

    ## add slopes/intercepts/etc.
    fits <- fitseg(x,y,segments)
    segments <- slopes(x,y,segments,fits)
    if ( "dpseg" %in% class(dpseg) ) dpseg$segments <- segments
    else dpseg <- segments
    dpseg
}

## get lm fits for final segments
## to be used also for generic predict
fitseg <- function(x,y,segments) {

    ## TODO: use try to catch and store official errors!
    fits <- rep(list("fiterror"), nrow(segments))

    for ( i in seq_len(nrow(segments)) ) {
        rng <- segments[i,"start"]:segments[i,"end"]
        yn <- y[rng]
        xn <- x[rng]
        ## remove NA or infinite data
        xn <- xn[!is.na(yn)&is.finite(yn)]
        yn <- yn[!is.na(yn)&is.finite(yn)]
        if ( length(yn)>1 ) {
            df <- data.frame(x=xn, y=yn)
            fits[[i]] <- lm(y ~ x, data=df) # defined herein!
        }
    }
    fits
}

## TODO: recognize and process fit types here!
slopes <- function(x, y, segments, fits)  {
    
    ## TODO: make this cleaner and fail-safer, recognize type of fit,
    ## and if known, get relevant parameters!
    ## TODO: handle fit errors via try!?
    error <- unlist(lapply(fits, function(x) x[[1]][1]=="fiterror"))

    ## get all fit coefficients
    res <- do.call(rbind,lapply(fits[!error], coef))
    
    ## recognize default linear fit?
    if (  "(Intercept)"==colnames(res)[1] )
        colnames(res) <- c("intercept.lm","slope.lm")
    
    resm <- matrix(NA, nrow=length(fits), ncol=ncol(res))
    colnames(resm) <- colnames(res)
    resm[!error,] <- res

    ## TODO: catch "essentially perfect fit" warning 
    
    ## get r.squared if provided in fit summary, otherwise just NULL
    if ( "r.squared" %in%  unlist(lapply(fits[!error], function(x)
        names(summary(x)))) ) {
        r2 <- matrix(NA, nrow=length(fits), ncol=1)
        r2[!error,] <- do.call(rbind,lapply(fits[!error],
                                            function(f) summary(f)$r.squared))
        resm <- cbind.data.frame(resm, r2.lm=r2)
    }

    ## get variance of residuals, if present
    var <- NULL
    if ( "residuals" %in% do.call(rbind,lapply(fits[!error], names)) ) {
        var <- matrix(NA, nrow=length(fits), ncol=1)
        var[!error,] <- do.call(rbind,lapply(fits[!error],
                                             function(f) var(f$residuals)))
        resm <- cbind.data.frame(resm, var.lm=var)
    }

    ## add coefficients and result summaries to data
    ## TODO: only cbind if present, fit-type specific
    segments <- cbind.data.frame(segments, resm)
    segments
}



#' Print method for linear segmentation result from \code{\link{dpseg}}.
#'
#' Prints the main result table \code{x$segments}, segment coordinates
#' and indices, and parameters from the recursion. See \code{\link{dpseg}}
#' for details.
#' @param x result object returned by function \code{\link{dpseg}}
#' @param ... further arguments to \code{print.data.frame}
#' @export
print.dpseg <- function(x, ...) {
    cat(paste("\ndynamic programming-based segmentation of",
             length(x$xy$x),"xy data points:\n\n"))
    print(x$segments, ...)
    parms <- paste(paste(names(x$parameters),
                         x$parameters,sep=": "),collapse="; ")
    cat(paste("\nParameters:", parms, "\n"))
}

#' Predict method for `dpseg` segmentations
#' 
#' Predicted values based on a data segmentation model from
#' \code{\link{dpseg}}.
#' @param object result object returned by function
#'     \code{\link{dpseg}}
#' @param xout new x-values at which to predict \eqn{\hat{y}}
#' @param ... not used
#' @return Returns predicted linear segments as x,y coordinates
#'     (\code{grDevices::xy.coords}) at \code{xout}.
#' @examples
#' x <- oddata$Time
#' y <- log(oddata$A5)
#' segs <- dpseg(x=x, y=y, P=0.0001)
#' 
#' ## predict method
#' plot(x=x, y=y, pch=19, cex=0.5)
#' lines(predict(segs), col=2, lwd=2)
#' @export
predict.dpseg <- function(object, xout, ...) {

    if ( missing(xout) )
        xout <- object$xy$x

    xnew <- NULL
    ynew <- NULL
    idxs <- NULL

     
    segs <- object$segments
    #fits <- object$fits

    ## calculated values
    slpcol <- "slope"
    icpcol <- "intercept"
    ## fall back on lm-derived?
    if ( !slpcol%in%colnames(segs) ) {
        idx <- grep("slope", colnames(segs))[1]
        pfix <- sub("slope","", colnames(segs)[idx])
        slpcol <- paste0("slope", pfix)
        icpcol <- paste0("intercept", pfix)
    }        
    if ( !slpcol%in%colnames(segs) ) 
        stop("No slope/intercept found: can't predict values for segments.",
             " Use `addLm(dpseg, x=x, y=y)` to add the original",
             " data and linear regression coefficients.")
    
    for ( i in seq_len(nrow(segs)) ) {
        rng <- c(segs[i,"x1"], segs[i,"x2"])
        ## TODO: account for jump?
        idx <- which(xout>=rng[1] & xout<=rng[2])
        xrng <- xout[idx]
        ## TODO: rm lm fits and use internal
        ## #
        if ( length(xrng)>0 )
            ##if ( use.fit )  ## keep this for now
            ##    yrng <- predict(fits[[i]], newdata= data.frame(x=xrng), ...)
            ##else
            yrng <- segs[i,slpcol]*xrng + segs[i,icpcol] 
        else 
            yrng <- rep(NA, length(xrng))
        idxs <- c(idxs,idx)                         
        ynew <- c(ynew,yrng) 
    }
    ## generate results for all xout
    yout <- rep(NA, length(xout))
    yout[idxs] <- ynew

    ## convert to xy plot coordinates
    grDevices::xy.coords(x=xout,y=yout)
}

#' Plot method for a \code{\link{dpseg}} segmentation model.
#' @param x result object returned by function \code{\link{dpseg}}
#' @param delog plot exp(y)
#' @param col optional color vector for segments
#' @param main plot title, dpseg parameters will be plotted if missing
#' @param ... arguments passed to default \code{\link{plot}} function
#' @return Silently returns the \code{x$segments} table , with color values
#' added if they were missing in the input.
#' @export 
plot.dpseg  <- function(x, delog=FALSE, col, main, ...) {

    segs <- x$segments
    xy <- x$xy

    if ( !"numeric" %in% class(xy$y) ) {
        stop("The results do not contain the original data and can not",
             " be plotted. Use `addLm(dpseg, x=x, y=y)` to add the original",
             " data and linear regression coefficients.")
    }
    
    if ( delog ) xy$y <- exp(xy$y)

    ## colors: generate if not present as argument or column
    if ( missing(col) ) col <- seq_len(nrow(segs))
    if ( !"col"%in%colnames(segs) )
        segs <- cbind(segs,col=col)

    cols <- rep(8, length(xy$x))
    for ( i in seq_len(nrow(segs)) ) # note: start overwrites previous end
        cols[segs[i,"start"]:segs[i,"end"]] <- segs[i,"col"]

    ## main text: parameters
    if ( missing(main) )
        main <- paste(paste(names(x$parameters),
                            x$parameters,sep=": "),collapse="; ")
    plot(xy, col=cols, main=main, ...)
    
    ## print predict result as lines
    ## TODO: when is the test required? for generic runs w/o add.lm?
    test <- try(pr <- predict.dpseg(x), silent=TRUE)
    if ( !"try-error" %in% class(test) ) {
        if ( delog ) pr$y <- exp(pr$y)
        lines(pr)
    } else warning(test[1])

    abline(v=xy$x[segs[seq_len(nrow(segs)),"end"]])
    abline(v=xy$x[segs[seq_len(nrow(segs)),"start"]],col=2,lty=2)
    #segs[["colors"]] <- cols
    invisible(segs)
}

#' Visualizes the \code{\link{dpseg}} segmentation recursion as a movie.
#'
#' Generates a movie of the calculation steps \eqn{j=1,\dots,n} while
#' looping through the recursion \eqn{S_j}. Plots are sent to the
#' active plot device or, if \code{path} is specified, to a video file
#' \code{<path>/<file.name>.<format>} via a system call to Image
#' Magick's \code{convert}.
#' Saving to a file likely only works on Linux systems with Image
#' Magick installed and \code{convert} available in the \code{$PATH}
#' environment variable. \code{format} are formats available for
#' \code{convert}, eg. \code{format="gif"} or \code{format="mpeg"}.
#' See the \code{vignette("dpseg")} for details on the plotted data.
#' @param dpseg result object of class \code{\link{dpseg}} returned by
#'     function \code{\link{dpseg}}
#' @param fix.ylim fix the y-axis of the \code{score} function
#' @param frames x range to show as movie frames
#' @param delay delay between frames in seconds, between x11 plot
#'     updates or as argument \code{-delay} to the system call to
#'     Image Magick's \code{convert}
#' @param repeat.last repeat list frame this many times
#' @param path path where both temporary jpeg files and the final
#'     movie file will be generated. If not specified the indidividual
#'     frames will be plotted to the active plot device.
#' @param file.name name of the generated video file
#'     <path>/<file.name>.<format>
#' @param format format of the video, all outputs that image magick's
#'     convert can generate, e.g. "mpg" or "gif"
#' @param res resolution of the generated movie (pixels per inch)
#' @param ylab left y-axis label, for the scoring funtion
#' @param ylab2 right y-axis label, for the original data
#' @param xlab x-axis label
#' @param ... arguments passed to default \code{\link{plot}} function
#' @examples
#' 
#' ## NOTE: requires that dpseg is run with store.matrix=TRUE
#' segs <- dpseg(x=oddata$Time, y=log(oddata$A3), minl=5, P=0.0001, store.matrix=TRUE)
#' 
#' ## View the algorithm in action:
#' movie(segs, delay=0)
#'
#' ## NOTE: if Image Magick's convert is installed you can set the path
#' ## option to save the movie as <path>/<file.name>.<format>, where format
#' ## can be "gif", "mpeg" or else, depending on the Image Magick installation.
#' @export
movie <- function(dpseg, fix.ylim=TRUE, frames, delay=0.1, repeat.last=5, 
                  ylab="scoring function", ylab2="y", xlab="x",
                  path, file.name="dpseg_movie", format="gif", res=200,...) {

    x <- dpseg$xy$x
    y <- dpseg$xy$y
    minl <- dpseg$parameters$minl
    maxl <- dpseg$parameters$maxl
    imax <- dpseg$imax
    imax[imax<1] <- 1 # TODO: fix upstream?
    segs <- dpseg$segments

    if ( missing(frames) )
        frames <- seq_along(x)

    have.matrix <- TRUE
    if ( is.null(dpseg$SCR) ) have.matrix <- FALSE
    else if (nrow(dpseg$SCR)==0 ) have.matrix <- FALSE
    if ( !have.matrix ) {
        warning("movie without scoring values; ",
                "run dpseg with store.matrix=TRUE for nicer movies")
        fix.ylim <- TRUE
        ylim <- 0:1
    } else {
        SCR <- dpseg$SCR
        SCR[is.infinite(SCR)] <- NA
        ylim <- range(SCR,na.rm=TRUE)
    }

    if ( !missing(path) )
        grDevices::jpeg(file.path(path,paste0(file.name,"_%04d.jpg")),
                        units="in",height=4,width=5,res=res)
    cat(paste("generating movie\n"))
    for ( j in c((minl+1):length(x),rep(length(x),repeat.last)) ) {

        if (!j%in%frames) next

        if ( !fix.ylim ) ylim <- range(SCR[,j],na.rm=TRUE)

        ## ensure to reset parameters
        opar <- par(no.readonly =TRUE)
        on.exit(par(opar))
        ## plot fitscore and segment borders
        par(mai=c(.5,.5,.3,.5),mgp=c(1.3,.2,0),tcl=-.25)
        plot(x,rep(0,length(x)), col=NA,
             ylab=NA, xlab=xlab, axes=FALSE, ylim=ylim, ...)
        axis(1)
        par(new=TRUE)
        if ( have.matrix ) {
            plot(x,SCR[,j], ylim=ylim, ylab=ylab, xlab=NA, axes=FALSE)
            axis(2)
        } 
        abline(v=x[j],lwd=3) # current position
        lcut <- j-c(maxl,minl) # minl/maxl cutoffs
        abline(v=x[lcut[lcut>0]],col=2)
        abline(v=x[imax[j]],col=4,lwd=3) # imax
        arrows(x0=x[imax[j]], x1=x[j], y0=ylim[1], lwd=2) # current segment
        par(new=TRUE)
        ## plot original data
        plot(x,y,axes=FALSE,col="gray",ylab=NA,xlab=NA)
        points(x[imax[j]:j], y[imax[j]:j],col=4)
        ## color already done segments
        sgs <- segs[segs[,"end"]<=imax[j],c("start","end"),drop=FALSE]
        if ( nrow(sgs)>=1 )
            for ( i in seq_len(nrow(sgs)) )
                points(x[sgs[i,1]:sgs[i,2]],y[sgs[i,1]:sgs[i,2]],col=i)
        axis(4,col="gray",col.axis="gray")
        mtext(ylab2,4, par("mgp")[1], col="gray")
        if ( j==length(x) ) ## final segment
            points(x[segs[nrow(segs),1]:segs[nrow(segs),2]],
                   y[segs[nrow(segs),1]:segs[nrow(segs),2]],col=nrow(segs))

        if( missing(path) ) Sys.sleep(delay)
    }

    if ( !missing(path) ) {
        dev.off()
        cmd <- paste0("convert -delay ", delay, " ",
                      file.path(path,paste0(file.name,"*.jpg")), " ",
                      file.path(path,paste0(file.name,".", format)))
        cat(paste("calling:\n", cmd, "\n"))
        system(cmd) # call Image Magick convert
        unlink(file.path(path,paste0(file.name,"*.jpg"))) # delete temp. files
        cat(paste("animation stored at:",
                  file.path(path,paste0(file.name,".", format)), "\n"))
    }
}

