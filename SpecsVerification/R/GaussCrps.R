#' Calculate the Continuous Ranked Probability Score (CRPS) for forecasts issued as Normal distributions
#'
#' @param mean A vector of length N. The forecast means.
#' @param sd A vector of length N. The forecast standard deviations.
#' @param obs A numeric vector of length N of real-valued verifying observations
#' @return numeric vector of length N with the CRPS values
#' @examples
#' data(eurotempforecast)
#' mean <- rowMeans(ens)
#' sd <- apply(ens, 1, sd)
#' mean(GaussCrps(mean, sd, obs))
#' @seealso EnsCrps, DressCrps, ScoreDiff, SkillScore
#' @references Gneiting et al (2005). Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics and Minimum CRPS Estimation. Mon. Wea. Rev. \doi{10.1175/MWR2904.1}
#' @export

GaussCrps <- function(mean, sd, obs) {

  # check if all vectors are of equal length
  stopifnot(length(mean) == length(obs))
  stopifnot(length(sd) == length(obs))

  # sample size
  N <- length(obs)

  # initialize crps vector, will be NA whenever sd < 0
  crps <- rep(NA, N)

  # for sd = 0, the normal is a step function, and the crps reduces to the
  # absolute difference
  i.zero <- which(sd == 0)
  if (length(i.zero) > 0) {
    crps[i.zero] <- abs(mean[i.zero] - obs[i.zero])
  }

  # crps for sd > 0
  i.pos <- which(sd > 0)
  if (length(i.pos) > 0) {
    z <- (obs - mean) / sd
    crps[i.pos] <- (sd * (z * (2 * pnorm(z) - 1) + 
                    2 * dnorm(z) - 1/sqrt(pi)))[i.pos]
  }

  # return crps vector
  return(crps)
}



