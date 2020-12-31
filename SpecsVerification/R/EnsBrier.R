#' Calculate the ensemble-adjusted Brier Score 
#'
#' @rdname EnsBrier
#' @param ens a N*R matrix representing N time instances of R-member ensemble forecasts of binary events; ens[t,r]=1 if the r-th ensemble member at time t predicted the event, otherwise ens[t,r]=0 
#' @param obs a numeric vector of length N with binary observations; obs[t]=1 if the event happens at time t, otherwise obs[t]=0 
#' @param R.new ensemble size for which the scores should be adjusted
#' @return numeric vector of length N with the ensemble-adjusted Brier scores 
#' @details `FairBrier(ens, obs)` returns `EnsBrier(ens, obs, R.new=Inf)`
#' @examples
#' data(eurotempforecast)
#' mean(EnsBrier(ens.bin, obs.bin, R.new=Inf))
#' @seealso EnsRps, EnsCrps, ScoreDiff, SkillScore
#' @references Ferro CAT, Richardson SR, Weigel AP (2008) On the effect of ensemble size on the discrete and continuous ranked probability scores. Meteorological Applications. \doi{10.1002/met.45}
#' @export

EnsBrier <- function(ens, obs, R.new=NA) {

  stopifnot(is.matrix(ens), 
            nrow(ens) == length(obs), 
            length(R.new)==1)

  # count number of ensemble members that predict the event
  i <- rowSums(ens==1)

  # calculate ensemble size as the number of valid forecasts, i.e. members that
  # are not NA|NaN|Inf
  R <- rowSums(is.finite(ens))

  # for R.new == NA (the default), no correction for ensemble size is performed
  if (length(R.new) == 1 & is.na(R.new[1])) {
    R.new <- R
  }

  # calculate ensemble-adjusted brier score
  br <- (i / R - obs)^2 - i * (R - i) / R / (R - 1) * (1 / R - 1 / R.new)
  
  # return the vector of brier scores
  return(br)

}



#' @rdname EnsBrier
#' @export
FairBrier <- function(ens, obs) {
  return(EnsBrier(ens, obs, R.new=Inf))
}

