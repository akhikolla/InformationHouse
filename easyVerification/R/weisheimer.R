# weisheimer.R Compute Reliability Categories as in Weisheimer et al. (2014)
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
#' Compute Reliability Categories as in Weisheimer et al. (2014)
#' 
#' This function implements the reliability categorization for forecasts of binary 
#' events as documented in Weisheimer et al. (2014). It has only been implemented for
#' category forecasts with categories defined relative to the forecast and observed
#' climatological distribution (i.e. without systematic bias). 
#' 
#' @param ens n x k matrix of n forecasts from k ensemble members
#' @param obs n verifying observations
#' @param pthresh probability threshold to convert to category forecasts. 
#'   If negative, event falling below threshold is used, else, event
#'   above threshold is used.
#' @param nboot number of bootstrap replicates to estimate 75 percent
#'   confidence interval
#' @param brier.thresholds Thresholds used to bin the forecasts (see 
#'   \code{\link[verification]{brier}})
#' @param ... additional arguments for compatibility with other scores
#' 
#' @keywords utilities
weisheimer <- function(ens, obs, pthresh=2/3, nboot=100, 
                       brier.thresholds=seq(0,1,0.2), ...) {
  if (requireNamespace("verification", quietly = TRUE)) {
    
    stopifnot(length(pthresh) == 1)
    prob <- abs(pthresh)
    ff <- convert2prob(ens, prob=prob)[,1.5 + 0.5*sign(pthresh)] / ncol(ens)
    oo <- convert2prob(obs, prob=prob)[,1.5 + 0.5*sign(pthresh)]
    stopifnot(length(ff) == length(oo))
    pp <- ifelse(pthresh > 0, 1 - prob, prob)
    
    bootfun <- function(ind=seq(along=ff)){
      ## compute brier score
      bs <- as.data.frame(verification::brier(oo[ind], ff[ind],
                                thresholds=brier.thresholds)[c('y.i', 'obar.i', 'prob.y')])
      ## normalise to go through origin
      bs$obar.i <- bs$obar.i - pp
      bs$y.i <- bs$y.i - pp
      
      ## compute weighted regression
      coef(lm(obar.i ~ y.i - 1, bs, weights=bs$prob.y))
    }
    
    ## set up output
    aa <- rep(NA, 3)
    aa[2] <- bootfun()
    
    ## run bootstrap
    aa.boot <- sapply(1:nboot, function(i){
      ind <- sample(seq(along=oo), length(oo), replace=TRUE)
      bootfun(ind)
    })
    
    ## compute lower and upper confidence bounds
    aa[c(1,3)] <- quantile(aa.boot, c(0.125, 0.875))
    
    ## compute categories
    out <- c(cat=1) + 
      as.numeric(aa[2] > 0) + 
      as.numeric(aa[1] > 0) + 
      as.numeric(aa[2] > 0.5 & aa[1] > 0) + 
      as.numeric(aa[3] > 1 & aa[1] > 0.5 & aa[1] < 1)
    
    return(out)
  } else {
    warning("package 'verification' must be installed for this function")
    return(NA)
  }
}
  
  