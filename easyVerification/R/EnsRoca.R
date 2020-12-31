# EnsRoca.R Area Under the ROC Curve
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

#' @name EnsRoca
#' 
#' @aliases oldEnsRoca EnsRocss
#'
#' @title Area Under the ROC Curve
#' 
#' @description Computes the area under the ROC curve given the observations. 
#'   \code{EnsRoca} computes the Area Under the Curve (AUC). For ease of interpretation,
#'   \code{EnsRocss} converts the AUC to the range from -1 to 1 with zero indicating
#'   a forecast with no discrimination.
#' 
#' @param ens n x j matrix of n probability forecasts for j categories
#' @param obs n x j matrix of occurence of n verifying observations in j categories
#' @param use.easy logical, should implementation of standard errors as implemented
#'   in \code{easyVerifcation} be used (see below)?
#' 
#' @section Standard Error: 
#'   
#'   If used with \code{SpecsVerification >= 0.5}, the standard errors as implemented
#'   in the function \code{SpecsVerification::Auc} are used. 
#'   
#'   If \code{use.easy = TRUE} or when used with an older version of \code{SpecsVerification}, 
#'   the standard error \eqn{\sigma} of the ROC area 
#'   skill score is given by the following formula after Broecker (2012).
#'   
#'   \deqn{\sigma^2 = \frac{1}{3} \left(\frac{1}{N_0} + \frac{1}{N_1} + 
#'   \frac{1}{N_0 N_1} \right)}{\sigma^2 = 1/3 (1/N0 + 1/N1 + 1/(N0 N1))}
#'   
#'   Where \eqn{\sigma} is the standard error, \eqn{N_1}{N1} the number of 
#'   events, and \eqn{N_0}{N0} the number of non-events in category \code{i}. 
#'   
#' @references Broecker, J. (2012). Probability forecasts. Forecast 
#'   Verification: A Practitioner's Guide in Atmospheric Science, Second 
#'   Edition, 119-139.
#' 
#' 
#' @examples
#' tm <- toymodel()
#' 
#' ## compute ROC area for tercile forecasts using veriApply
#' veriApply("EnsRoca", fcst=tm$fcst, obs=tm$obs, prob=1:2/3)
#' 
#' @seealso \code{\link{veriApply}}, \code{\link{EnsRocss}}
#' 
#' @export
EnsRoca <- function(ens, obs, use.easy=FALSE){
  stopifnot(is.matrix(ens), is.matrix(obs), length(obs) == length(ens))
  if (packageVersion("SpecsVerification") >= 0.5 & ! use.easy){
    ens.prob <- count2prob(ens, type=4)
    aucfun <- get("Auc")
    roc.area <- lapply(1:ncol(ens.prob), function(i) {
      auc <- try(aucfun(ens.prob[,i], obs[,i]), silent=TRUE)
      if (class(auc) == 'try-error') auc <- as.list(as.numeric(rep(NA, 2)))
      return(auc)})
    roc.area <- as.list(unlist(roc.area, recursive=FALSE))
    names(roc.area) <- paste0(rep(paste0('cat', 1:ncol(ens)), each=2), c("", ".sigma"))
    ## resort for consistency with old approach
    roc.area <- roc.area[c(seq(1, length(roc.area), 2), seq(2, length(roc.area), 2))]
  } else {
    roc.area <- EnsRocaCpp(ens, obs)
    roc.area <- as.list(roc.area)
    names(roc.area) <- paste0('cat', seq(along=roc.area))
    ## compute sigma
    N1 <- apply(obs, 2, sum)
    N0 <- nrow(obs) - N1
    roc.sigma <- sqrt(1/12 * (1/N0 + 1/N1 + 1/(N0*N1)))
    roc.sigma[N1 == 0] <- NA
    roc.sigma <- as.list(roc.sigma)
    names(roc.sigma) <- paste0('cat', seq(along=roc.sigma), '.sigma')
    roc.area <- c(roc.area, roc.sigma)    
  }
  return(roc.area)
}

#' @rdname EnsRoca
#' @export
EnsRocss <- function(ens, obs, use.easy=FALSE){
  roc.area <- EnsRoca(ens=ens, obs=obs, use.easy=use.easy)
  roc.skill <- lapply(roc.area, function(x) 2*x)
  notsigma <- -grep("sigma", names(roc.skill))
  roc.skill[notsigma] <- lapply(roc.skill[notsigma], function(x) x - 1)
  return(roc.skill)
}


oldEnsRoca <- function(ens, obs){
  stopifnot(is.matrix(ens), is.matrix(obs), length(obs) == length(ens))
  rs <- rowSums(ens)
  if (any(rs != 1)){
    ## convert number of occurences to probabilities
    ens <- ens / rs
  }
  n.event <- colSums(obs)
  n.total <- nrow(obs)
  ens.rank <- apply(ens, 2, rank)
  mean.rank <- colSums(ens.rank * obs) / n.event
  roc.area <- (mean.rank - (n.event + 1)/2) / (n.total - n.event)
  roc.area[n.event == 0] <- NA
  roc.area <- as.list(roc.area)
  names(roc.area) <- paste0('cat', seq(along=roc.area))
  return(roc.area)
}