# EnsCorr.R Correlation with Ensemble Mean
#
#     Copyright (C) 2016 MeteoSwiss
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' Calculate Fair Ranked Probability Skill Score Against Climatological
#' Reference Forecast.
#' 
#' Calculate the fair ranked probability skill score (fair RPSS) between an
#' ensemble forecasts and a climatological reference forecast derived from the
#' observations. The categories of the climatological reference forecast have
#' been defined based on the distribution of the observations and the exact
#' forecast probabilities are known. The 'fair' correction therefore should not
#' be applied to the reference forecast.
#' 
#' @param ens	N*K matrix. ens[i,j] is the number of ensemble members that
#'   predict category j at time i.
#' @param ens.ref N*K matrix, similar to ens
#' @param obs N*K matrix. obs[i,j] = 1 if category j is observed at time i, 0
#'   otherwise.
#' @param format additional argument for use with \code{SpecsVerification >= 0.5}.
#'   Do not change this argument manually (except when using \code{climFairRpss}, 
#'   as standalone function).
#'   
#' @return A list with the following elements: \code{rpss|skillscore}: The value of the
#' skill score. \code{sigma.rpss|skillscore.sd}: The standard deviation of the skill score,
#' approximated by propagation of uncertainty. Please note that the naming changes with the
#' new version of \code{SpecsVerification}.
#'   
#' @examples
#' tm <- toymodel()
#' 
#' ## compute RPSS using veriApply
#' veriApply("climFairRpss", tm$fcst, tm$obs, prob=1:2/3)
#' 
#' @seealso \code{\link{veriApply}}
#'   
#' @export
climFairRpss <- function (ens, ens.ref, obs, format=c("category", "member")) {
  if (packageVersion("SpecsVerification") >= 0.5){
    skillfun <- 'SkillScore'
    out <- get(skillfun)(changearg(EnsRps, format=format)(ens, obs, R.new=Inf),
                         changearg(EnsRps, format=format)(ens.ref, obs, R.new=NA),
                         handle.na="use.pairwise.complete")    
  } else {
    stopifnot(all(dim(ens) == dim(obs)))
    stopifnot(all(dim(ens.ref) == dim(obs)))
    rps.ens <- FairRps(ens, obs)
    rps.ref <- EnsRps(ens.ref, obs)
    i.nna <- which(!is.na(rps.ens + rps.ref))
    rps.ens <- rps.ens[i.nna]
    rps.ref <- rps.ref[i.nna]
    N <- length(i.nna)
    rpss <- 1 - mean(rps.ens)/mean(rps.ref)
    rpss.sigma <- 1/sqrt(N) * sqrt(var(rps.ens)/mean(rps.ref)^2 + 
                                     var(rps.ref) * mean(rps.ens)^2/mean(rps.ref)^4 - 2 * 
                                     cov(rps.ens, rps.ref) * mean(rps.ens)/mean(rps.ref)^3)
    out <- list(rpss = rpss, rpss.sigma = rpss.sigma)
    
  }
  
  return(out)
}