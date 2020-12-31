#' Creating Time-Series according to example model
#' 
#' Creates a time series according to the block model with independent innovations.
#' @param rand.gen optional: a function to generate the innovations
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a vector consisting of a realization of the signal model
#' @details See Appendix B in the reference for details about the signal model.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal
modelSignal.blocks <- function(rand.gen=rnorm, ...) {
  signal.blocks.means <- c(0, 14.64, -3.66, 7.32, -7.32,
                          10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)
  signal.blocks.Cpts <- c(c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)-1, 2048)
  signal.blocks.lengths <- c(204, diff(signal.blocks.Cpts))
  signal.blocks.sigma <- 10
  mu_t <- numeric(0)
  for(i in 1:length(signal.blocks.lengths)) {
    mu_t <- c(mu_t, rep(signal.blocks.means[i],
                        signal.blocks.lengths[i]))
  }
  sigma_t <- rep(signal.blocks.sigma, length(mu_t))
  return(list(mu_t=mu_t,
              sigma_t=sigma_t))
}

#' Creating Time-Series according to example model
#' 
#' Creates a time series according to the block model with independent innovations.
#' @param rand.gen optional: a function to generate the innovations
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a vector consisting of a realization of the signal model
#' @details See Appendix B in the reference for details about the signal model.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal
modelSignal.fms <- function(rand.gen=rnorm, ...) {
  signal.fms.Cpts<-c(c(139,226,243,300,309,333)-1,497)
  signal.fms.sigma <- 0.3
  signal.fms.means<-c(-0.18,0.08,1.07,-0.53,0.16,-0.69,-0.16)
  signal.fms.lengths <- c(signal.fms.Cpts[1], diff(signal.fms.Cpts))
  mu_t <- numeric(0)
  for(i in 1:length(signal.fms.lengths)) {
    mu_t <- c(mu_t, rep(signal.fms.means[i],
                        signal.fms.lengths[i]))
  }
  sigma_t <- rep(signal.fms.sigma, length(mu_t))
  return(list(mu_t=mu_t,
              sigma_t=sigma_t))
}

#' Creating Time-Series according to example model
#' 
#' Creates a time series according to the block model with independent innovations.
#' @param rand.gen optional: a function to generate the innovations
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a vector consisting of a realization of the signal model
#' @details See Appendix B in the reference for details about the signal model.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal
modelSignal.mix <- function(rand.gen=rnorm, ...) {
  signal.mix.Cpts<-c(10,20,40,60,90,120,160,200,250,300,360,420,490,560)
  signal.mix.sigma <- 4
  signal.mix.means<-c(7,-7,6,-6,5,-5,4,-4,3,-3,2,-2,1,-1)
  signal.mix.lengths <- c(10, diff(signal.mix.Cpts))
  mu_t <- numeric(0)
  for(i in 1:length(signal.mix.lengths)) {
    mu_t <- c(mu_t, rep(signal.mix.means[i],
                        signal.mix.lengths[i]))
  }
  sigma_t <- rep(signal.mix.sigma, length(mu_t))
  return(list(mu_t=mu_t,
              sigma_t=sigma_t))
}

#' Creating Time-Series according to example model
#' 
#' Creates a time series according to the block model with independent innovations.
#' @param rand.gen optional: a function to generate the innovations
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a vector consisting of a realization of the signal model
#' @details See Appendix B in the reference for details about the signal model.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal
modelSignal.teeth10 <- function(rand.gen=rnorm, ...) {
  signal.teeth.Cpts<-c(10,20,30,40,50,60,70,80,90,100,110,120,130,140)
  signal.teeth.sigma <- 0.4
  signal.teeth.means<-c(0,1,0,1,0,1,0,1,0,1,0,1,0,1)
  signal.teeth.lengths <- c(signal.teeth.Cpts[1], diff(signal.teeth.Cpts))
  mu_t <- numeric(0)
  for(i in 1:length(signal.teeth.lengths)) {
    mu_t <- c(mu_t, rep(signal.teeth.means[i],
                        signal.teeth.lengths[i]))
  }
  sigma_t <- rep(signal.teeth.sigma, length(mu_t))
  return(list(mu_t=mu_t,
              sigma_t=sigma_t))
}

#' Creating Time-Series according to example model
#' 
#' Creates a time series according to the block model with independent innovations.
#' @param rand.gen optional: a function to generate the innovations
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a vector consisting of a realization of the signal model
#' @details See Appendix B in the reference for details about the signal model.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal
modelSignal.stairs10 <- function(rand.gen=rnorm, ...) {
  signal.stairs.Cpts<-c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
  signal.stairs.sigma <- 0.3
  signal.stairs.means<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  signal.stairs.lengths <- c(signal.stairs.Cpts[1], diff(signal.stairs.Cpts))
  mu_t <- numeric(0)
  for(i in 1:length(signal.stairs.lengths)) {
    mu_t <- c(mu_t, rep(signal.stairs.means[i],
                        signal.stairs.lengths[i]))
  }
  sigma_t <- rep(signal.stairs.sigma, length(mu_t))
  return(list(mu_t=mu_t,
              sigma_t=sigma_t))
}