#' Test data with piecewise constant mean
#' 
#' Generate piecewise stationary time series with independent innovations and change-points in the mean.
#' @param model a string indicating from which model a realisation is to be generated;
#' possible values are "custom" (for user-specified model
#' using \code{lengths}, \code{means} and \code{sds}), and
#' "blocks", "fms", "mix", "stairs10", "teeth10" (for the referenced test signals)
#' @param lengths use iff \code{model = "custom"}; an integer vector for the lengths of the piecewise stationary segments
#' @param means use iff \code{model = "custom"}; a numeric vector for the means of the piecewise stationary segments
#' @param sds use iff \code{model = "custom"}; a numeric vector for the deviation scaling of the piecewise stationary segments.
#' The values are multiplied to the outcome of \code{rand.gen}, coinciding with the standard
#' deviation in the case of standard normal innovations (\code{rand.gen = rnorm})
#' @param rand.gen optional; a function to generate the noise/innovations
#' @param seed optional; if a seed value is provided (\code{!is.null(seed)}), 
#' then \code{set.seed(seed)} is called beforehand)
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return a list containing the following entries:
#' \itemize{
#' \item{x}{ a numeric vector containing a realisation of the piecewise time series model, 
#' given as signal + noise}
#' \item{mu}{ mean vector of piecewise stationary time series model }
#' \item{sigma}{ scaling vector of piecewise stationary time series model }
#' \item{cpts}{ a vector of change-points in the piecewise stationary time series model }
#' }
#' @details See Appendix B in the reference for details about the test signals.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @examples
#' # visualise estimated changepoints by solid vertical lines
#' # and true changepoints by broken vertical lines
#' td <- testData(lengths = c(50, 50, 200, 300, 300), means = c(0, 1, 2, 3, 2.3), 
#' sds = rep(1, 5), seed = 123)
#' mbu <- multiscale.bottomUp(td$x)
#' plot(mbu, display = "data")
#' abline(v = td$cpts, col = 2, lwd = 2, lty = 2)
#' 
#' # visualise estimated piecewise constant signal by solid line
#' # and true signal by broken line
#' td <- testData("blocks", seed = 123)
#' mlp <- multiscale.localPrune(td$x)
#' plot(mlp, display = "data")
#' lines(td$mu, col = 2, lwd = 2, lty = 2)
#' @importFrom stats rnorm
#' @export
testData <- function(model=c('custom', 'blocks', 'fms', 'mix', 'stairs10', 'teeth10')[1],
                                           lengths=NULL, means=NULL, sds=NULL,
                                           rand.gen=rnorm, seed=NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  signal <- testSignal(model=model, lengths=lengths, means=means, sds=sds)
  n <- length(signal$mu_t)
  ts <- signal$mu_t + rand.gen(n, ...)*signal$sigma_t
  return(list(x = ts, mu = signal$mu_t, sigma = signal$sigma_t, cpts = which(diff(signal$mu_t) != 0)))
}

#' Piecewise constant test signal
#' 
#' Produce vectors of mean and dispersion values for generating piecewise stationary time series.
#' @param model a string indicating from which model a realisation is to be generated;
#' possible values are "custom" (for user-specified model
#' using \code{lengths}, \code{means} and \code{sds}), and
#' "blocks", "fms", "mix", "stairs10", "teeth10" (for the referenced test signals)
#' @param lengths use iff \code{model = "custom"}; an integer vector for the lengths of the piecewise stationary segments
#' @param means use iff \code{model = "custom"}; a numeric vector for the means of the piecewise stationary segments
#' @param sds use iff \code{model = "custom"}; a numeric vector for the deviation scaling of the piecewise stationary segments.
#' @return a list containing the following entries:
#' \itemize{
#'   \item{mu_t}{ mean vector of piecewise stationary model time series}
#'   \item{sigma_t}{ deviation scaling vector of piecewise stationary model time series}
#' }
#' @details See Appendix B in the reference for details about the test signals.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @keywords internal 
testSignal <- function(model=c('custom', 'blocks', 'fms', 'mix', 'stairs10', 'teeth10')[1],
                                       lengths=NULL, means=NULL, sds=NULL) {
  if (model=='blocks') {
    res <- modelSignal.blocks()
  } else if (model=='fms') {
    res <- modelSignal.fms()
  } else if (model=='mix') {
    res <- modelSignal.mix()
  } else if (model=='stairs10') {
    res <- modelSignal.stairs10()
  } else if (model=='teeth10') {
    res <- modelSignal.teeth10()
  } else if (model=='custom') {
    stopifnot(!is.null(lengths))
    stopifnot(!is.null(means))
    stopifnot(!is.null(sds))
    stopifnot(length(lengths)==length(means))
    stopifnot(length(lengths)==length(sds))
    mu_t <- numeric(0)
    sigma_t <- numeric(0)
    for(i in 1:length(lengths)) {
      mu_t <- c(mu_t, rep(means[i], lengths[i]))
      sigma_t <- c(sigma_t, rep(sds[i], lengths[i]))
    }
    res <- list(mu_t=mu_t,
                sigma_t=sigma_t)
  } else {
    stop('Unkonwn model string')
  }
  return(res)
}
