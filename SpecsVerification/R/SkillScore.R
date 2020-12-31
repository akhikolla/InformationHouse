#' Calculate a skill score and assess uncertainty. 
#'
#' A skill score is defined as (mean score - mean reference score) / (perfect score - mean reference score). The skill score is zero if the mean score of the forecast equals the mean score of the reference forecast, and equals one if the mean score of the forecast equals the best possible score. Uncertainty is assessed by estimating the standard deviation of the skill score by propagation of uncertainty.
#'
#' @param scores vector of verification scores
#' @param scores.ref vector of verification scores of the reference forecast, must be of the same length as `scores`
#' @param N.eff user-defined effective sample size to be used to estimate the sampling uncertainty; if NA, the length of `scores` is used; default: NA
#' @param score.perf a numeric constant, indicating the value that the score would assign to the perfect forecast
#' @param handle.na how should missing values in scores vectors be handled; possible values are 'na.fail' and 'use.pairwise.complete'; default: 'na.fail'
#' @return vector with skill score and its estimated standard deviation
#' @examples
#' data(eurotempforecast)
#' SkillScore(EnsCrps(ens, obs), EnsCrps(ens[, 1:2], obs))
#' @seealso ScoreDiff
#' @export
SkillScore <- function(scores, scores.ref, N.eff=NA, score.perf=0, handle.na=c('na.fail', 'use.pairwise.complete')) {

  ## sanity checks
  N.eff <- N.eff[1L]
  score.perf <- score.perf[1L]
  stopifnot(N.eff > 0 | is.na(N.eff))
  stopifnot(length(scores) == length(scores.ref))
  stopifnot(is.finite(score.perf))


  ## handle NA's
  handle.na = match.arg(handle.na)
  if (handle.na == "na.fail") {
    if (any(is.na(c(scores, scores.ref)))) {
      stop("missing values")
    }
  } 
  if (handle.na == "use.pairwise.complete") {
    nna <- !is.na(scores) & !is.na(scores.ref)
    if (all(nna == FALSE)) {
      stop("there are no complete pairs of scores")
    }
    scores <- scores[nna]
    scores.ref <- scores.ref[nna]
  }


  ## _after_ removing any missing values, deal with user-defined effective sample size
  if (is.na(N.eff)) {
    N.eff <- length(scores)
  }




  ## calculations

  # calculate mean scores, shift by score.perf
  score <- mean(scores) - score.perf
  score.ref <- mean(scores.ref) - score.perf

  # calculate skill score
  skillscore <- 1 - score / score.ref

  # calculate auxiliary quantities
  v.score <- var(scores)
  v.score.ref <- var(scores.ref)
  cov.score   <- cov(scores, scores.ref)

  # calculate skill score standard deviation by error propagation 
  sqrt.na <- function(z) {
    z[z<0] <- NA
    return(sqrt(z))
  }
  skillscore.sigma <- 
         1 / sqrt(N.eff) * sqrt.na( v.score / score.ref^2 + 
         v.score.ref * score^2 / score.ref^4 - 
         2 * cov.score * score / score.ref^3)

  # set skillscore.sigma to NA if not finite
  if (!is.finite(skillscore.sigma)) {
    skillscore.sigma <- NA 
  }

  #return
  c(skillscore=skillscore, skillscore.sd=skillscore.sigma)

}





##################################################################
# the following functions are included for backward compatibility with
# SpecsVerification version < 1.0.0



#' Calculate DressCrps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param dressed.ens the ensemble
#' @param dressed.ens.ref the reference ensemble
#' @param obs the observation
#' @return DressCrps Skill Score
#' @seealso SkillScore DressCrps DressEnsemble
#' @export

DressCrpss <- function(dressed.ens, dressed.ens.ref, obs) {
  SkillScore(DressCrps(dressed.ens, obs), DressCrps(dressed.ens.ref, obs), handle.na="use.pairwise.complete")
}



#' Calculate EnsBrier Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param tau not used
#' @return EnsBrier skill score
#' @seealso SkillScore EnsBrier
#' @export

EnsBrierSs <- function(ens, ens.ref, obs, tau=NA) {
  SkillScore(EnsBrier(ens, obs), EnsBrier(ens.ref, obs), handle.na="use.pairwise.complete")
}





#' Calculate EnsCrps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @return EnsCrps skill score
#' @seealso SkillScore EnsCrps
#' @export

EnsCrpss <- function(ens, ens.ref, obs) {
  SkillScore(EnsCrps(ens, obs), EnsCrps(ens.ref, obs), handle.na="use.pairwise.complete")
}







#' Calculate EnsRps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param format see `EnsRps`
#' @return EnsRps skill score
#' @seealso SkillScore EnsRps
#' @export

EnsRpss <- function(ens, ens.ref, obs, format=c('category', 'members')) {
  SkillScore(EnsRps(ens=ens, obs=obs, format=format), 
             EnsRps(ens=ens.ref, obs=obs, format=format), 
             handle.na="use.pairwise.complete")
}





#' Calculate FairBrier Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param tau not used
#' @return FairBrier skill score
#' @seealso SkillScore EnsBrier
#' @export

FairBrierSs <- function(ens, ens.ref, obs, tau=NA) {
  SkillScore(EnsBrier(ens, obs, R.new=Inf), EnsBrier(ens.ref, obs, R.new=Inf), handle.na="use.pairwise.complete")
}





#' Calculate FairCrps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @return FairCrps skill score
#' @seealso SkillScore EnsCrps
#' @export

FairCrpss <- function(ens, ens.ref, obs) {
  SkillScore(EnsCrps(ens, obs, R.new=Inf), EnsCrps(ens.ref, obs, R.new=Inf), handle.na="use.pairwise.complete")
}







#' Calculate FairRps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param ens the ensemble 
#' @param ens.ref the reference ensemble
#' @param obs the observation
#' @param format see `EnsRps`
#' @return FairRps skill score
#' @seealso SkillScore EnsRps
#' @export

FairRpss <- function(ens, ens.ref, obs, format=c('category', 'members')) {
  SkillScore(EnsRps(ens, obs, R.new=Inf, format=format), 
             EnsRps(ens.ref, obs, R.new=Inf, format=format), 
             handle.na="use.pairwise.complete")
}




#' Calculate GaussCrps Skill Score (deprecated, use function SkillScore instead)
#'
#' @param mean forecast means
#' @param sd forecast standard deviations
#' @param mean.ref reference forecast means
#' @param sd.ref reference forecast standard deviations
#' @param obs the observation
#' @return GaussCrps skill score
#' @seealso SkillScore GaussCrps
#' @export

GaussCrpss <- function(mean, sd, mean.ref, sd.ref, obs) {
  SkillScore(GaussCrps(mean, sd, obs), GaussCrps(mean.ref, sd.ref, obs), handle.na="use.pairwise.complete")
}




