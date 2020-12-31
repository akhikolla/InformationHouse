#' Calculate the squared error between forecast and observation
#'
#' @param fcst a N-vector representing N time instances of real-valued forecasts
#' @param obs a N-vector representing N time instances of real-valued observations
#' @return numeric N-vector of squared errors
#' @examples
#' data(eurotempforecast)
#' mean(SqErr(rowMeans(ens), obs))
#' @seealso AbsErr, ScoreDiff, SkillScore
#' @export

SqErr <- function(fcst, obs) {

  # calculate squared difference between forecasts and observations
  sq.err <- (fcst - obs) ^ 2

  # return 
  return(sq.err)

}


