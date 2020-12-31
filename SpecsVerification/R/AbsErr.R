#' Calculate the absolute error between forecast and observation
#'
#' @param fcst a N-vector representing N time instances of real-valued forecasts
#' @param obs a N-vector representing N time instances of real-valued observations
#' @return numeric N-vector of absolute errors |fcst - obs|
#' @examples
#' data(eurotempforecast)
#' mean(AbsErr(rowMeans(ens), obs))
#' @seealso SqErr, ScoreDiff, SkillScore
#' @export

AbsErr <- function(fcst, obs) {

  # calculate absolute difference between forecasts and observations
  abs.err <- abs(fcst - obs)

  # return 
  return(abs.err)

}

