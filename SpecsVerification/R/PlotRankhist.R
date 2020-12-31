#' Plotting function for rank histograms
#' @description Plots a rank histogram in different modes.
#' @param rank.hist A vector or rank counts.
#' @param mode Either "raw" (default) or "prob.paper". Whether to draw the raw rank histogram, or the rank histogram on probability paper. See Details.
#' @details The plotting modes currently implemented are:
#'
#' raw (the default): A simple bar plot of the counts provided by the `rank.hist` argument.
#'
#' prob.paper: The individual counts given by `rank.hist` are transformed to their cumulative probabilities under the binomial distribution with parameters `N` and `1/K`, where `N=sum(rank.hist)` and `K=length(rank.hist)`. This transformation makes possible an assessment of the observed rank counts under the hypothesis of equally likely ranks. The y-axis on the left indicates the cumulative probabilities. The intervals on the right of the plot indicate central 90, 95, and 99 percent _simultaneous_ confidence intervals. That is, if all ranks were equally likely on average, approximately 90 percent of all rank histograms would be _completely_ contained in the 90 percent interval and approximately 10 percent of all rank histograms would have _at least_ one bar that falls outside this interval.
#' @examples
#' data(eurotempforecast)
#' rank.hist <- Rankhist(ens, obs)
#' PlotRankhist(rank.hist, mode="prob.paper")
#' @references Anderson J.L. (1996). A Method for Producing and Evaluating Probabilistic Forecasts from Ensemble Model Integrations. J. Climate, 9, 1518--1530.
#' Broecker J. (2008). On reliability analysis of multi-categorical forecasts. Nonlin. Processes Geophys., 15, 661-673.
#' @seealso Rankhist, TestRankhist
#' @export

PlotRankhist <- function(rank.hist, mode="raw") {
#
#
# Plot the rank histogram raw or on probability paper
#
# Usage: PlotRankhist(rh, mode="raw")
#
# Arguments:
#   * rank.hist ... a vector of rank counts (see function `rankhist()`)
#   * mode      ... whether to plot raw or on probability paper to 
#                   assess flatness
#
# Details:
#
#   * the `raw` mode simply plots a barplot of the rank histogram counts
#   * the `prob.paper` mode transforms the observed rank histogram counts to
#     cumulative probabilities under Binomial(N, 1/(K+1)), plots them on a logit
#     scale, and adds simultaneous consistency intervals
#
# Return value:
#   * none (yet)
#
# Author: 
#
#    Stefan Siegert 
#    s.siegert@exeter.ac.uk 
#    October 2013
#
# Example:
# 
#   ens <- matrix(rnorm(500), 100, 5)
#   obs <- rnorm(100)
#   rh <- rankhist(ens, obs)
#   PlotRankhist(rank.hist = rh, mode="prob.paper")
#
# References: 
#   Broecker 2008 http://dx.doi.org/10.5194/npg-15-661-2008
#
################################
  if (mode == "prob.paper") {
    N <- sum(rank.hist)
    K <- length(rank.hist) - 1
    # cumulative binomial likelihood of the observed rank counts
    nuh <- pbinom(q=rank.hist, size=N, prob=1/(K+1))
    # log-odds ratio
    lornuh <- log(nuh / (1 - nuh))
    # clip very small and large bars
    clipval <- 8
    i.clip.min <- which(lornuh < -clipval)
    lornuh[i.clip.min] <- -clipval
    i.clip.max <- which(lornuh >  clipval)
    lornuh[i.clip.max] <-  clipval
    # prepare for plotting
    offs <- 8
    bar.wd <- 0.9
    par(oma=c(0, 0, 0, 4), cex.lab=0.8, cex.axis=0.8)
    plot(NULL, xlim=c(0,K+2), ylim=c(0,2*offs), axes=F, xlab="rank i", ylab="")
    b <- barplot(lornuh + offs, add=T, axes=FALSE, width=bar.wd, col=gray(0.5))
    points(b[i.clip.min], lornuh[i.clip.min]+offs, pch=25, bg="black", cex=.6)
    points(b[i.clip.max], lornuh[i.clip.max]+offs, pch=24, bg="black", cex=.6)
    # axis labels
    xlabels <- paste(1:(K+1))
    axis(1, at=b, labels=xlabels)
    pvals <- c(0.001,0.01,0.1,0.5,0.9,0.99,0.999)
    yvals <- log(pvals/(1-pvals))
    axis(2, at=yvals+offs, labels=pvals, las=1)
    # confidence intervals, corrected for multiple testing
    ci.str <- c("90%","95%","99%")
    k <- 0
    for (p in c(0.9,0.95,0.99)) {
      k <- k+1
      del.b <- diff(b)[1]
      ci <- c(1-p^(1/(K+1)),p^(1/(K+1)))
      ci <- log(ci/(1-ci))+offs
      axis(4,at=ci,labels=c("",""),tcl=0.5,line=k)
      axis(4,at=log(c(0.1,0.9)/(1-c(0.1,0.9)))+offs,labels=c("",""),col="white",lwd=3,line=k)
      mtext(ci.str[k],side=4,padj=-1,line=k,cex=0.8)
      lines(x=c(b[1]-.7*bar.wd,b[K+1]+0.7*bar.wd),y=rep(ci[1],2),lty=2,xpd=TRUE)
      lines(x=c(b[1]-.7*bar.wd,b[K+1]+0.7*bar.wd),y=rep(ci[2],2),lty=2,xpd=TRUE)
    }
  } else if (mode == "raw") {
    par(mar=c(3, 3, 1, 0), cex.lab=0.8, cex.axis=0.8)
    bp <- barplot(rank.hist, xlab=NA, ylab=NA, col=gray(0.5), axes=FALSE)
    axis(1, at=bp, line=.2, labels=paste(1:length(rank.hist)))
    axis(2, las=1)
    mtext(side=1, text="rank i", line=2, cex=.8)
    mtext(side=2, text="count", line=2, cex=.8)
  }
}


