### EXPERIMENTAL: parameter estimation and scans

## estimate good range for P, when using variance of residuals,
## NOTE: taking 1/100 of this worked for OD4 data from 20200629_topA
#' Estimate a starting value for penalty P.
#'
#' The break-point penalty P in a \code{\link{dpseg}} recursion,
#' should be in the range of expected values of the scoring
#' function. To find a good initial estimate for P when using the
#' default scoring fuction (see \code{\link{dpseg}}), the data is
#' smoothed by \code{smooth.spline} and the variance of residuals
#' reported.
#' @param x x-values
#' @param y y-values
#' @param plot plot the \code{\link{smooth.spline}}
#' @param ... parameters for \code{\link{smooth.spline}}
#' @return Returns a double, variance of residuals of a spline fit
#'     (\code{var(smooth.spline(x,y, ...)$y -y)})
#' @examples
#' x <- oddata$Time
#' y <- log(oddata$A5)
#' p <- estimateP(x=x, y=y, plot=TRUE)
#' plot(dpseg(x=x, y=y, jumps=TRUE, P=round(p,3)))
#' @export
estimateP <- function(x, y, plot=FALSE, ...) {
    rm <- is.na(y)|is.infinite(y)
    x <- x[!rm]
    y <- y[!rm]
    ys <- smooth.spline(x,y, ...)$y
    if ( plot ) {
        plot(x, y, cex=1, pch=19, col="#00000099")
        lines(x, ys, col="red", lwd=1.5)
    }
    var(ys -y)
}

#' Scan over different penalty \code{P} values
#'
#' Runs the \code{\link{dpseg}} recursion for different values of the
#' penalty parameter \code{P} and returns a matrix with the used
#' \code{P} values, the resulting number of segments and (optionally)
#' the median of segment variance of residuals.
#' @param x x-values
#' @param y y-values
#' @param P vector of penalties P to scan
#' @param var add the median of the variances of residuals of all
#'     segments to output (save time by \code{var=FALSE})
#' @param use.matrix use the stored scoring function matrix for more
#'     efficient scans; set this to \code{FALSE} if you run into
#'     memory problems
#' @param plot plot results
#' @param verb print progress messages
#' @param ... parameters for \code{\link{dpseg}} (except \code{P})
#' @return Returns a matrix with the penalties P in the first column,
#'     the number of segments in the second column and the median of
#'     variances in the third column.
#' @examples
#' x <- oddata$Time
#' y <- log(oddata$A5)
#' par(mai=c(par("mai")[1:3], par("mai")[2])) # to show right axis
#' sp <- scanP(x=x, y=y, P=seq(-.01,.1,length.out=50), plot=TRUE)
#'@export
scanP <- function(x, y, P, var=TRUE,
                  use.matrix=TRUE, plot=TRUE, verb=1, ...) {

    ps <- P
    nsg <- rep(NA, length(ps))

    if ( var )
        vars <- nsg
    
    if ( verb>0 )
        cat(paste("running dpseg", length(ps), "times: "))
    
    ## if enough memory is available, we can run this
    ## more efficiently by storing the scoring matrix
    if ( use.matrix ) {
        ## run once with store-matrix
        dp <- dpseg(x=x, y=y, P=ps[1], store.matrix=TRUE, verb=0, ...)

        if ( verb>0 ) cat(".")

        nsg[1] <- nrow(dp$segments)
        if ( var )
            vars[1] <- median(dp$segments$var)
        
        ## and rerun using pre-calculated matrix
        if ( length(ps)>1 )
            for ( i in 2:length(ps) ) {
                dpi <- dpseg(y=dp$SCR, P=ps[i], store.matrix=FALSE, verb=0, ...)
                nsg[i] <- nrow(dpi$segments)
                if ( var ) {
                    dpi <- addLm(dpi, x=x, y=y)
                    vars[i] <-  median(dpi$segments$var.lm)
                }
                if ( verb>0 ) cat(".")
            }
    } else {
        for ( i in seq_along(ps) ) {
            dpi <- dpseg(x=x, y=y, P=ps[i],
                         store.matrix=FALSE, verb=0, ...)
            nsg[i] <- nrow(dpi$segments)
            if ( var )
                vars[i] <-  median(dpi$segments$var)
            if ( verb>0 ) cat(".")
        }
    }
    if ( verb>0 ) cat("\n")
    res <- cbind(P=ps,"number of segments"=nsg)
    if ( var )
        res <- cbind(res, "median of variances"=vars)
    if ( plot ) {
        ## ensure to reset parameters
        opar <- par(no.readonly =TRUE)
        on.exit(par(opar))
        
        plot(res, type="l")
        points(res, pch=19, cex=.5)
        if ( var ) {
            par(new=TRUE)
            plot(res[,1],res[,3], type="l", col=2,axes=FALSE,xlab=NA,ylab=NA)
            points(res[,1], res[,3], col=2, pch=19, cex=.5)
            axis(4,col=2,col.axis=2)
            mtext("median of variances", 4, par("mgp")[1], col=2)
        }
            
    }
    invisible(res)
}

