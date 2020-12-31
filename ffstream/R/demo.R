#' Demo for ffstream
#'
#' Provides a demonstration of the AFF method detecting changepoints in a 
#' stream.
#'
#' @param showPlot Boolean flag; if \code{TRUE}, then a plot is generated.
#'                 Default is \code{FALSE}.
#'
#' @param returnStream Boolean flag; if \code{TRUE}, then return the stream 
#'                     as part of the list returned by the demo.
#'                     Default is \code{FALSE}.
#'
#' @param plotSmall Boolean flag; if \code{TRUE}, creates a small plot, as
#'                  needed for the vignette. Default is \code{FALSE}.
#'
#' @details This method generates a stream with three changepoints, and finds 
#'          the changepoints with AFF. Also creates a plot of the data and 
#'          the changepoints if the \code{showPlot} flag is set to \code{TRUE}.
#'          The observations are shown in black, the true changepoints are
#'          shown as red dotted vertical lines, and the detected (estimated)
#'          changepoints are shown as blue dashed lines. The following is 
#'          returned in a list:
#'          \describe{
#'                    \item{tau}{The location of the true changepoints.}
#'
#'                    \item{tauhat}{The detected (estimated) changepoints.}
#'
#'                    \item{method}{The method used, in this case AFF.}
#'
#'                    \item{param}{The data frame with the parameters used in 
#'                                 the AFF method, in this case, 
#'                          \describe{
#'                                    \item{alpha}{The significance level,}
#'
#'                                    \item{eta}{The step size in the 
#'                                               gradient descent, whose
#'                                               value is not particularly
#'                                               important,}
#'
#'                                    \item{BL}{The length of the burn-in
#'                                              period.}
#'
#'
#'                                   }
#'                                 }
#'                   }
#'
#'
#'
#' @section Author:
#' Dean Bodenham
#'
#'
#' @section References:
#' D. A. Bodenham and N. M. Adams (2016) 
#' \emph{Continuous monitoring for changepoints in data 
#' streams using adaptive estimation}. 
#' Statistics and Computing  
#' doi:10.1007/s11222-016-9684-8
#'
#'
#' @examples
#' df <- demo_ffstream()
#'
#' \dontrun{
#' demo_ffstream(showPlot=TRUE)
#'}
#'
#'
#' @export
demo_ffstream <- function(showPlot=FALSE, returnStream=FALSE, plotSmall=FALSE){

    seednum <- 3
    regimeLength <- 150
    numChanges <- 3
    tau <- seq_len(numChanges) * regimeLength
    stream <- makeStreamMeanChangeR(seednum=seednum, 
                                    regimeLength=regimeLength, 
                                    numChanges=numChanges)
    eta <- 0.001
    BL <- 50
    alpha <- 0.01
    AFFList <- detectAFFMean(x=stream, alpha=alpha, eta=eta, BL=BL)
    if (showPlot){

        #part of moving legend to bottom
        #bottom left top right
#        par(oma = c(0, 0, 0, 0))


        xlab <- "Observation \n"
        ylab <- "Value"
        maintitle <- "ffstream demo"
        tauCol <- "red"
        tauLty <- 3
        tauhatLty <- 2
        tauhatCol <- "blue"
        lineSize <- 2
        figSizeWidth <- 4
        figSizeHeight <- 3
        legend1 <- "True changepoints (tau)"
        legend2 <- "Detected changepoints (tauhat)"
        legendfontsize <- 0.8
        legendx <- -100
        legendy <- 7.5
        fontsize=1

        if (plotSmall){
            figSizeWidth <- 2
            figSizeHeight <- 1.5
            legendfontsize <- 0.5
            legendx <- -200
            legendy <- 8.5
            lineSize <- 1.5
        }

        par(pin=c(figSizeWidth, figSizeHeight))
        plot(x=seq_len(length(stream)), y=stream, type='l', col="black",
                xlab=xlab, ylab=ylab, main=maintitle, cex=fontsize)
        abline(v=tau, col=tauCol, lty=tauLty, lwd=lineSize)
        abline(v=AFFList$tauhat, col=tauhatCol, lty=tauhatLty, lwd=lineSize)

        legendvec <- c(legend1, legend2)
        colours <- c(tauCol, tauhatCol)
        legendlty <- c(tauLty, tauhatLty)
        legend(legendx, legendy, legendvec, xpd = TRUE, horiz = TRUE, 
                inset = c(0, 0), bty = "n", lty=legendlty, 
                col = colours, lwd = lineSize, cex=legendfontsize)
    } 

    param <- data.frame(alpha=alpha, eta=eta, BL=BL)
    if (returnStream){
        returnList <- list(stream=stream, tau=tau, tauhat=AFFList$tauhat, 
                           method="AFF", param=param)
    } else {
        returnList <- list(tau=tau, tauhat=AFFList$tauhat, 
                           method="AFF", param=param)
    }
    return(returnList)
}
