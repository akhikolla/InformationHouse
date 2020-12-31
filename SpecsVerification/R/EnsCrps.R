#' Calculate the ensemble-adjusted Continuous Ranked Probability Score (CRPS)
#'
#' @rdname EnsCrps
#' @param ens a N*R matrix representing N time instances of real-valued R-member ensemble forecasts
#' @param obs a numeric vector of length N with real-valued observations
#' @param R.new positive number, can be `Inf`, ensemble size for which the scores should be adjusted, default is NA for no adjustment
#' @return numeric vector of length N with the ensemble-adjusted CRPS values
#' @details `FairCrps(ens, obs)` returns `EnsCrps(ens, obs, R.new=Inf)`
#' @examples
#' data(eurotempforecast)
#' mean(EnsCrps(ens, obs, R.new=Inf))
#' @seealso EnsBrier, EnsRps, DressCrps, GaussCrps, ScoreDiff, SkillScore
#' @references Ferro CAT, Richardson SR, Weigel AP (2008) On the effect of ensemble size on the discrete and continuous ranked probability scores. Meteorological Applications. \doi{10.1002/met.45}
#' @export
EnsCrps = function(ens, obs, R.new=NA) {

  stopifnot(is.matrix(ens), 
            nrow(ens) == length(obs), 
            length(R.new)==1)

  crps = enscrps_cpp(ens, obs, R.new)

  # return the vector of crps
  return(crps)

}


#' @rdname EnsCrps
#' @export
FairCrps = function(ens, obs) {
  return(EnsCrps(ens, obs, R.new=Inf))
}
