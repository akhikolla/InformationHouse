#' Calculate correlation difference between a forecast and a reference forecast, and assess uncertainty
#'
#' @param fcst vector of forecasts
#' @param fcst.ref vector of reference forecasts
#' @param obs vector of observations
#' @param N.eff user-defined effective sample size to be used in hypothesis test and for confidence bounds; if NA, the length of `obs` is used after removing missing values; default: NA
#' @param conf.level confidence level for the confidence interval; default = 0.95
#' @param handle.na how should missing values in forecasts and observations be handled; possible values are 'na.fail' and 'only.complete.triplets'; default: 'na.fail'
#' @return vector with correlation difference, one-sided p-value, and central confidence interval at the user-defined confidence level
#' @examples
#' data(eurotempforecast)
#' CorrDiff(rowMeans(ens), ens[, 1], obs)
#' @seealso Corr
#' @references Steiger (1980): Tests for comparing elements of a correlation matrix. Psychological Bulletin. \doi{10.1037/0033-2909.87.2.245} 
#' Zou (2007): Toward using confidence intervals to compare correlations. Psychological Methods. \doi{10.1037/1082-989X.12.4.399}
#' @export
CorrDiff <- function(fcst, fcst.ref, obs, N.eff=NA, conf.level=0.95, handle.na="na.fail") {

  ## sanity checks
  stopifnot(length(fcst) == length(obs))
  stopifnot(length(fcst.ref) == length(obs))


  ## handle NA's
  if (handle.na == "na.fail") {
    if (any(is.na(c(fcst, fcst.ref, obs)))) {
      stop("missing values")
    }
  } else if (handle.na == "only.complete.triplets") {
    nna <- !is.na(fcst) & !is.na(fcst.ref) & !is.na(obs)
    if (all(nna == FALSE)) {
      stop("there are no complete sets of forecasts and observations")
    }
    fcst <- fcst[nna]
    fcst.ref <- fcst.ref[nna]
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
  

  ## calculate correlation coefficients and their confidence intervals
  cc.fcst <- cor(fcst, obs)
  cc.ref <- cor(fcst.ref, obs)


  ## calculate correlation difference
  cc.diff <- cc.fcst - cc.ref


  ## auxiliary quantities
  r12 <- cc.fcst
  r13 <- cc.ref
  r23 <- cor(fcst, fcst.ref)


  ## confidence interval, according to zou 2007, example 2, fail if N <= 3
  if (conf.level <= 0 | conf.level >= 1) {
    conf.level <- NA
  }
  if (is.na(conf.level) | N <= 3) {
    L <- U <- NA
  } else {
    # individual confidence limits of cc.ens and cc.ref
    z.fcst <- 0.5 * log((1 + cc.fcst) / (1 - cc.fcst))
    z.ref <- 0.5 * log((1 + cc.ref) / (1 - cc.ref))

    alpha <- (1. - conf.level) / 2.
    l1 <- tanh(atanh(r12) + qnorm(alpha)/sqrt(N-3))
    u1 <- tanh(atanh(r12) + qnorm(1-alpha)/sqrt(N-3))
    l2 <- tanh(atanh(r13) + qnorm(alpha)/sqrt(N-3))
    u2 <- tanh(atanh(r13) + qnorm(1-alpha)/sqrt(N-3))

    # correlation between the two corcoefs r12 and r13
    c.12.13 <- 
      ((r23 - 0.5 * r12 * r13) * (1 - r12*r12 - 
      r13*r13 - r23*r23) + r23*r23*r23) / ((1 - r12*r12) * 
      (1 - r13*r13))
    # lower confidence limit
    L <- r12 - r13 - sqrt((r12-l1)^2 + (u2 - r13)^2 - 
         2*c.12.13*(r12-l1)*(u2-r13))
    # upper confidence limit
    U <- r12 - r13 + sqrt((u1 - r12)^2 + (r13-l2)^2 - 
         2*c.12.13*(u1-r12)*(r13-l2))
  }


  ## p value of one-sided test for equality of dependent correlation
  ## coefficients (steiger 1980 Eq 7)
  if (N <= 3) {
    p.value <- NA
  } else {
    R <- (1-r12*r12-r13*r13-r23*r23) + 2*r12*r13*r23
    t <- (r12 - r13) * sqrt((N-1)*(1+r23) / (2 * ((N-1)/(N-3)) * 
      R + 0.25*(r12+r13)^2 * (1-r23)^3))
    p.value <- 1 - pt(t, df=N-3)
  }


  ## return
  ret <- c(cc.diff, p.value, L, U)
  names(ret) <- c("corr.diff", "p.value", "L", "U")
  return(ret)
}

