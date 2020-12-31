#' Calculate the ensemble-adjusted Quadratic Score (QS) for categorical forecasts
#'
#' @rdname EnsQs
#' @param ens a N*R matrix of integers, representing N time instances of categorical ensemble forecasts; ens[t,r] indicates the category index that the r-th ensemble member forecasts at time t
#' @param obs a vector of length N, obs[t] is the category that occurred at time t
#' @param R.new ensemble size for which the scores should be adjusted
#' @return numeric vector of length N with the ensemble-adjusted quadratic score values
#' @details `FairQs(ens, obs)` returns `EnsQs(ens, obs, R.new=Inf)`
#' @examples
#' data(eurotempforecast)
#' EnsQs(ens.cat, obs.cat, R.new=Inf)
#' @details It is assumed that the smallest class index is 1, and the largest class index is calculated by max(c(ens,obs))
#' @seealso EnsBrier, EnsRps, EnsCrps, ScoreDiff, SkillScore
#' @export
EnsQs <- function(ens, obs, R.new=NA) {


  ## calculate histogram at each time, the tabulate function automatically removes NAs and NaNs
  K <- max(c(ens,obs), na.rm=TRUE)
  ens.hist <- t(apply(ens, 1, tabulate, nbins=K))
  obs.hist <- t(sapply(obs, tabulate, nbins=K))


  ## calculate time length, and ensemble sizes at each time
  N <- nrow(ens.hist)
  R <- rowSums(ens.hist)


  ## calculate unadjusted quadratic score, make use of row-first indexing when
  ## dividing the matrix ens.hist by the vector R
  q.score <- rowSums((ens.hist / R - obs.hist) ^ 2)


  ## check if ensemble-adjustment is requested by user
  R.new <- R.new[1]
  if (!is.na(R.new)) {
    if (R.new < 2) {
      q.score <- q.score + NA
    } else {
      adjustment <- rowSums(-1 * (1/R - 1/R.new) * ens.hist * (R - ens.hist) / R / (R-1))
      # catch one-member ensembles, which cannot be adjusted
      adjustment[!is.finite(adjustment)] <- NA
      q.score <- q.score + adjustment
    }
  }

  return(q.score)
}

#' @rdname EnsQs
#' @export
FairQs = function(ens, obs) {
  return(EnsQs(ens, obs, R.new=Inf))
}
