#' Reliability diagram for probability forecasts
#'
#' @param probs vector of N probability forecasts for the event obs=1
#' @param obs vector of N binary observations, event/no event are coded as 0/1
#' @param bins binning to estimate the calibration function (see Details), default: 10
#' @param nboot number of bootstrap resamples to calculate the consistency bars, default: 500
#' @param plot logical, whether to plot the reliability diagram, default: FALSE
#' @param plot.refin Whether to add the frequency distribution of the forecasts to the reliability diagram. default: TRUE
#' @param cons.probs The width of the consitency intervals. default: 0.95. Set to NA for no consistency bars.
#' @param attributes locical, whether attributes lines are included in the diagram. default: FALSE
#' @param handle.na how should missing values be handled; possible values are 'na.fail' and 'use.pairwise.complete'; default: 'na.fail'
#' @return a data.frame with nrows equal to the number of bins (given by the `bins` argument), with columns: average forecast probability per bin, conditional event frequency per bin, lower and upper limit of the consistency bar per bin, number of forecast probabilities per bin, lower and upper bin limit
#' @examples
#' data(eurotempforecast)
#' p <- rowMeans(ens.bin)
#' ReliabilityDiagram(p, obs.bin, plot=TRUE)
#' @details
#' To estimate the reliability curve, the unit line is categorised into discrete bins, provided by the `bins` argument. If `bins` is a single number, it specifies the number of equidistant bins. If `bins` is a vector of values between zero and one, these values are used as the bin-breaks. 
#' @references 
#' Jolliffe IT, Stephenson DB, eds. (2012): Forecast verification: A practitioner's guide in atmospheric science. John Wiley & Sons, 2012. ISBN: 978-0-470-66071-3
#' Broecker J, Smith LA (2007): Increasing the Reliability of Reliability Diagrams. Wea. Forecasting, 22, 651--661 \doi{10.1175/WAF993.1}
#' @export

ReliabilityDiagram <- 
function(probs, obs, bins=10, nboot=500, 
         plot=FALSE, plot.refin=TRUE, 
         cons.probs=0.95, attributes=FALSE, 
         handle.na=c("na.fail", "use.pairwise.complete"))
{

  handle.na = match.arg(handle.na)
  if (handle.na == "na.fail") {
    if (any(is.na(c(probs, obs)))) {
      stop("Missing values in probs and/or obs. Try setting handle.na=\"use.pairwise.complete\"")
    }
  }

  # sanity checks
  stopifnot(length(probs) == length(obs))
  stopifnot(nboot >= 0)
  stopifnot(is.logical(plot), is.logical(plot.refin))
  stopifnot(all(obs %in% c(0,1,NA)))

  # for backward compatibility
  if (length(cons.probs) > 1) {
    cons.probs <- max(cons.probs) - min(cons.probs)
  }

  # filter complete forecast-observation pairs
  if (handle.na == "use.pairwise.complete") {
    nna <- !is.na(probs) & !is.na(obs)
    if (all(nna == FALSE)) {
      stop("there are no complete forecast-observation pairs")
    }
    probs <- probs[nna]
    obs <- obs[nna]
  }
  stopifnot(all(probs >= 0), all(probs <= 1))

  # some definitions and corrections
  n <- length(obs)
  nboot <- floor(nboot)


  #############################################
  # reliability analysis
  #############################################
  # estimate refinement function
  if (length(bins) == 1) {
    nbins <- floor(bins)
    brx <- seq(0, 1, length.out=nbins+1) 
  } else {
    nbins <- length(bins) - 1
    bins <- sort(bins)
    stopifnot(min(bins)<= 0 & max(bins) >= 1)
    brx <- bins
  }
  h <- hist(probs, breaks=brx, plot=FALSE)$counts        


  # estimate calibration function
  g <- hist(probs[obs==1], breaks=brx, plot=FALSE)$counts
  obar.i <- g / h 
  obar.i[ is.nan(obar.i) ] <- NA
  

  # calculate in-bin averages
  p.bins <- as.numeric(cut(probs, breaks=brx, include.lowest=TRUE))
  p.avgs <- sapply(seq(nbins), 
                   function(ii) mean(probs[p.bins == ii], na.rm=TRUE))
  p.avgs[ is.nan(p.avgs) ] <- NA



  #############################################
  # consistency resampling (broecker and smith 2007)
  #############################################
  if (nboot & !is.na(cons.probs)) {
    resamp.mat <- matrix(nrow=0, ncol=nbins)
    # the resampling function
    sample.rel.diag <- function() {
      p.hat <- sample(x=probs, size=n, replace=TRUE)
      x.hat <- rbinom(n=n, size=1, prob=p.hat)
      hh <- hist(p.hat, breaks=brx, plot=FALSE)$counts        
      gg <- hist(p.hat[x.hat==1], breaks=brx, plot=FALSE)$counts
      return(gg / hh)
    }
    l <- replicate(nboot, sample.rel.diag())
    resamp.mat <- t(l)
    alpha <- 1 - cons.probs
    cons.bars <- apply(resamp.mat, 2, 
                       function(z) quantile(z, c(alpha/2, 1-alpha/2), na.rm=TRUE))
  } else {
    cons.bars <- matrix(NA, ncol=nbins, nrow=2)
  }


  #############################################
  # plot the reliability diagram
  #############################################
  if (plot) {
    # reliability plot
    old.par <- par(no.readonly = TRUE) 
    on.exit(par(old.par))
    par(las=1)
    plot(NULL, xlim = c(0,1), ylim = c(0,1),
       xlab= "Forecast probability",
       ylab="Observed relative frequency", xaxs="i", yaxs="i")
   if (attributes) {
       obs.clim <- mean(obs)
       a   <- (1+obs.clim)/2 
       b   <- obs.clim / 2
       x.p <- c(obs.clim, obs.clim, 1, 1, 0, 0)
       y.p <- c(0, 1, 1, a, b, 0)
       polygon(x.p, y.p, col = "#e6e6e6", border=NA)
    }
    # consistency bars
    for (i in 1:length(p.avgs)) {
        lines(rep(p.avgs[i], 2), cons.bars[, i], col="#CCCCCC", lwd=6, lend=1)
    }
    # reliability points, diagonal, and box
    points(p.avgs, obar.i, col = "black", pch = 1, lwd=2, type="b")
    lines(c(0,1), c(0,1), lty=1)
    box()
    if (plot.refin) {
      # refinement histogram in lower corner
      pp<- par("plt")
      par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) , las=1)
      par(new = TRUE)
      barplot(h, axes = FALSE, axisnames = FALSE)
      axis(4)
      box() 
    }
  }

  # for return, calculate the lower and upper bin limits
  bin.lower <- brx[-length(brx)]
  bin.upper <- brx[-1]

  #############################################
  # return data
  #############################################
  ret.df <- data.frame(p.avgs=p.avgs, cond.probs=obar.i, 
                       cbar.lo=cons.bars[1,], cbar.hi=cons.bars[2,], 
                       p.counts=h, bin.lower=bin.lower, bin.upper=bin.upper)
  return(ret.df)
}


