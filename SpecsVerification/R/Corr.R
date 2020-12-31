#' Calculate correlation between forecasts and observations, and assess uncertainty
#'
#' @param fcst vector of forecasts
#' @param obs vector of observations
#' @param N.eff user-defined effective sample size to be used in hypothesis test and for confidence bounds; if NA, the length of `obs` is used after removing missing values; default: NA
#' @param conf.level confidence level used the confidence interval; default = 0.95
#' @param handle.na how should missing values in forecasts and observations be handled; possible values are 'na.fail' and 'use.pairwise.complete'; default: 'na.fail'
#' @return vector with correlation, one-sided p-value, and central confidence interval at the user-defined confidence level
#' @examples
#' data(eurotempforecast)
#' Corr(rowMeans(ens), obs)
#' @seealso CorrDiff
#' @references Von Storch, Zwiers (2001): Statistical analysis in climate research. Cambridge University Press.
#' @export
Corr <- function(fcst, obs, N.eff=NA, conf.level=0.95, handle.na="na.fail") {

  ## sanity checks
  stopifnot(length(fcst) == length(obs))


  ## handle NA's
  if (handle.na == "na.fail") {
    if (any(is.na(c(fcst, obs)))) {
      stop("missing values")
    }
  } else if (handle.na == "use.pairwise.complete") {
    nna <- !is.na(fcst) & !is.na(obs)
    if (all(nna == FALSE)) {
      stop("there are no complete sets of forecasts and observations")
    }
    fcst <- fcst[nna]
    obs <- obs[nna]
  } else {
    stop("unknown 'handle.na' argument")
  }


  ## define sample size; the case of N <= 3 etc is handeled by the testing functions
  N.eff <- N.eff[1L]
  if (!is.na(N.eff)) {
    N <- N.eff
  } else {
    N <- length(obs)
  }
  

  ## calculate correlation
  cc <- cor(obs, fcst)


  # calculate confidence interval by Fisher z-transform, return NA if N < 4
  alpha <- (1. - conf.level) / 2.
  probs <- c(alpha, 1 - alpha)
  if (N <= 3) {
    ci <- c(NA_real_, NA_real_)
  } else {
    ci <- tanh(atanh(cc) + qnorm(probs)/sqrt(N-3))
  }

  # calculate p value by one-sided t-test (vonstorch & zwiers, p149)
  if (N <= 2) {
    p.val <- NA
  } else {
    t <- sqrt((N-2) * cc * cc / (1 - cc * cc))
    p.val <- pt(t, df=N-2, lower.tail=FALSE)
  }
  ret <- c(cc, p.val, ci)
  names(ret) <- c("corr", "p.value", "L", "U")

  return(ret)
}

