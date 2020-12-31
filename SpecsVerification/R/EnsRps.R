#' Calculate the ensemble-adjusted Ranked Probability Score (RPS) for categorical forecasts
#'
#' @rdname EnsRps
#' @param ens matrix with N rows representing N time instances of categorical ensemble forecasts as follows: If `format = category` (the default), then ens[t,r] indicates the category that the r-th ensemble member predicts for time t. Note that categories must be positive integers. If `format = members`, then ens[t,k] is the number of ensemble members that predict category k at time t.
#' @param obs vector of length N, or matrix with N rows, representing the N observed category as follows: If `format = category', obs is a vector and obs[t] is the category observed at time t. If `format = members`, obs is a matrix where obs[t,k] = 1 (and zero otherwise) if category k was observed at time t
#' @param R.new ensemble size for which the scores should be adjusted, defaults to NA (no adjustment)
#' @param format string, 'category' (default) or 'members' (can be abbreviated). See descriptions of arguments `ens` and `obs` for details.
#' @return numeric vector of length N with the ensemble-adjusted RPS values
#' @details `FairRps(ens, obs)` returns `EnsRps(ens, obs, R.new=Inf)`
#' @examples
#' data(eurotempforecast)
#' EnsRps(ens.cat, obs.cat, R.new=Inf)
#' @seealso EnsBrier, EnsQs, EnsCrps
#' @export
EnsRps <- function(ens, obs, R.new=NA, format=c('category', 'members')) {

  format = match.arg(format)

  if (format == 'category' & any(c(ens, obs) <= 0)) {
    stop('Categories must be positive numbers.')
  }
  
  if (format == 'category') {
    ## calculate histogram at each time, the tabulate function automatically removes NAs and NaNs
    K <- max(c(ens,obs), na.rm=TRUE)
    ens.hist <- t(apply(ens, 1, tabulate, nbins=K))
    obs.hist <- t(sapply(obs, tabulate, nbins=K))
  } else {
    K <- ncol(ens)
    ens.hist <- ens
    obs.hist <- obs
  }

  ## calculate time length, and ensemble sizes at each time
  N <- nrow(ens.hist)
  R <- rowSums(ens.hist)

  # accumulate ensemble members and observations
  ens.cum <- t(apply(ens.hist, 1, cumsum))
  obs.cum <- t(apply(obs.hist, 1, cumsum))

  ## calculate unadjusted quadratic score, make use of row-first indexing when
  ## dividing the matrix ens.cum by the vector R
  rps <- rowSums((ens.cum / R - obs.cum) ^ 2)


  # check if ensemble-adjustment is requested by user
  R.new <- R.new[1]
  if (!is.na(R.new)) {
    if (R.new < 2) {
      rps <- rps + NA
    } else {
      adjustment <- rowSums(-1 * (1/R - 1/R.new) * ens.cum * (R - ens.cum) / R / (R-1))
      # catch one-member ensembles, which cannot be adjusted
      adjustment[!is.finite(adjustment)] <- NA
      rps <- rps + adjustment
    }
  }

  return(rps)
  
}

#' @rdname EnsRps
#' @export
FairRps <- function(ens, obs, format=c('category', 'members')) {
  return(EnsRps(ens, obs, R.new=Inf, format))
}
