#' Calculate average score difference and assess uncertainty
#'
#' Calculate the difference (mean score of the reference forecast) minus (mean score of the forecast). Uncertainty is assessed by the Diebold-Mariano test for equality of predictive accuracy.
#'
#' @param scores vector of verification scores
#' @param scores.ref vector of verification scores of the reference forecast, must be of the same length as `scores`
#' @param N.eff user-defined effective sample size to be used in hypothesis test and for confidence bounds; if NA, the length of `scores` is used; default: NA
#' @param conf.level confidence level for the confidence interval; default = 0.95
#' @param handle.na how should missing values in scores vectors be handled; possible values are 'na.fail' and 'use.pairwise.complete'; default: 'na.fail'
#' @return vector with mean score difference, estimated standard error of the mean, one-sided p-value of the Diebold-Mariano test, and the user-specified confidence interval
#' @examples
#' data(eurotempforecast)
#' ScoreDiff(EnsCrps(ens, obs), EnsCrps(ens[, 1:2], obs))
#' @seealso SkillScore
#' @references Diebold, Mariano (1995): Comparing Predictive Accuracy. Journal of Business & Economic Statistics. \url{https://www.jstor.org/stable/1392185}
#' @export
ScoreDiff <- function(scores, scores.ref, N.eff=NA, conf.level=0.95, handle.na="na.fail") {


  ## sanity checks
  N.eff <- N.eff[1L]
  stopifnot(N.eff > 0 | is.na(N.eff))
  stopifnot(length(scores) == length(scores.ref))


  ## handle NA's
  if (handle.na == "na.fail") {
    if (any(is.na(c(scores, scores.ref)))) {
      stop("missing values")
    }
  } else if (handle.na == "use.pairwise.complete") {
    nna <- !is.na(scores) & !is.na(scores.ref)
    if (all(nna == FALSE)) {
      stop("there are no complete pairs of scores")
    }
    scores <- scores[nna]
    scores.ref <- scores.ref[nna]
  } else {
    stop("unknown 'handle.na' argument")
  }


  ## _after_ removing any missing values, deal with user-defined effective sample size
  if (is.na(N.eff)) {
    N.eff <- length(scores)
  }


  ## calculations

  # calculate loss-differentials
  d <- scores.ref - scores

  # calculate mean score difference between reference forecast and forecast
  d.bar <- mean(d)

  # calculate variance of loss-differentials
  var.d <- var(d)
  
  # calculate standard error of the mean loss-differential, use effective sample size
  d.bar.sd <- sqrt(var.d / N.eff)

  # calculate p-value using asymptotic test by Diebold-Mariano (1995), return
  # NA if variance is non-positive
  if (d.bar.sd > 0) {
    p.value <- pnorm(d.bar / d.bar.sd, lower.tail=FALSE)
  } else {
    p.value <- NA
  }

  # calculate confidence interval of the mean
  if (conf.level <= 0 | conf.level >= 1) {
    conf.level <- NA
  }
  if (is.na(conf.level)) {
    ci <- c(NA, NA)
  } else {
    alpha <- (1. - conf.level) / 2.
    ci <- qnorm(p=c(alpha, 1-alpha), mean=d.bar, sd=d.bar.sd)
  }
  

  ## return vector including the score difference, error of the mean, p.value,
  # and confidence interval
  ret <- c(score.diff=d.bar, score.diff.sd=d.bar.sd, p.value=p.value, L=ci[1], U=ci[2])
  return(ret)

}



##################################################################
# the following functions are included for backward compatibility with
# SpecsVerification version < 1.0.0



#' Calculate DressCrps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param dressed.ens the ensemble
#' @param dressed.ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @return mean DressCrps difference
#' @seealso ScoreDiff DressCrps DressEnsemble
#' @export

DressCrpsDiff <- function(dressed.ens, dressed.ens.ref, obs, probs=NA) {
  ScoreDiff(DressCrps(dressed.ens, obs), DressCrps(dressed.ens.ref, obs), handle.na="use.pairwise.complete")
}


#' Calculate DressIgn Difference (deprecated, use function ScoreDiff instead)
#'
#' @param dressed.ens the ensemble
#' @param dressed.ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @return mean DressIgn difference
#' @seealso ScoreDiff DressIgn
#' @export

DressIgnDiff <- function(dressed.ens, dressed.ens.ref, obs, probs=NA) {
  ScoreDiff(DressIgn(dressed.ens, obs), DressIgn(dressed.ens.ref, obs), handle.na="use.pairwise.complete")
}




#' Calculate EnsBrier Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param tau not used
#' @param probs not used
#' @return mean EnsBrier difference
#' @seealso ScoreDiff EnsBrier
#' @export

EnsBrierDiff <- function(ens, ens.ref, obs, tau=NA, probs=NA) {
  ScoreDiff(EnsBrier(ens, obs), EnsBrier(ens.ref, obs), handle.na="use.pairwise.complete")
}





#' Calculate EnsCrps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @return mean EnsCrps difference
#' @seealso ScoreDiff EnsCrps
#' @export

EnsCrpsDiff <- function(ens, ens.ref, obs, probs=NA) {
  ScoreDiff(EnsCrps(ens, obs), EnsCrps(ens.ref, obs), handle.na="use.pairwise.complete")
}







#' Calculate EnsRps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @param format see `EnsRps`
#' @return mean EnsRps difference
#' @seealso ScoreDiff EnsRps
#' @export

EnsRpsDiff <- function(ens, ens.ref, obs, probs=NA, format=c("category", "members")) {
  ScoreDiff(EnsRps(ens=ens, obs=obs, R.new=NA, format=format), 
            EnsRps(ens=ens.ref, obs=obs, R.new=NA, format=format), 
            handle.na="use.pairwise.complete")
}





#' Calculate FairBrier Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param tau not used
#' @param probs not used
#' @return mean FairBrier difference
#' @seealso ScoreDiff EnsBrier
#' @export

FairBrierDiff <- function(ens, ens.ref, obs, tau=NA, probs=NA) {
  ScoreDiff(EnsBrier(ens, obs, R.new=Inf), EnsBrier(ens.ref, obs, R.new=Inf), handle.na="use.pairwise.complete")
}





#' Calculate FairCrps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @return mean FairCrps difference
#' @seealso ScoreDiff EnsCrps
#' @export

FairCrpsDiff <- function(ens, ens.ref, obs, probs=NA) {
  ScoreDiff(EnsCrps(ens, obs, R.new=Inf), EnsCrps(ens.ref, obs, R.new=Inf), handle.na="use.pairwise.complete")
}







#' Calculate FairRps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param probs not used
#' @param format see `EnsRps`
#' @return mean FairRps difference
#' @seealso ScoreDiff EnsRps
#' @export

FairRpsDiff <- function(ens, ens.ref, obs, probs=NA, format=c("category", "members")) {
  ScoreDiff(EnsRps(ens=ens, obs=obs, R.new=Inf, format=format), 
            EnsRps(ens=ens.ref, obs=obs, R.new=Inf, format=format), 
            handle.na="use.pairwise.complete")
}




#' Calculate GaussCrps Difference (deprecated, use function ScoreDiff instead)
#'
#' @param mean forecast means
#' @param sd forecast standard deviations
#' @param mean.ref reference forecast means
#' @param sd.ref reference forecast standard deviations
#' @param obs the observation
#' @param probs not used
#' @return mean GaussCrps difference
#' @seealso ScoreDiff GaussCrps
#' @export

GaussCrpsDiff <- function(mean, sd, mean.ref, sd.ref, obs, probs=NA) {
  ScoreDiff(GaussCrps(mean, sd, obs), GaussCrps(mean.ref, sd.ref, obs), handle.na="use.pairwise.complete")
}




