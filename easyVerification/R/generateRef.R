# generateRef Generate Probabilistic Climatological Ensemble Forecast from Observations
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

#' @name generateRef
#' @aliases indRef
#'   
#' @title Generate Probabilistic Climatological Ensemble Forecast from 
#'   Observations
#'   
#' @description To generate reference ensemble forecasts for forecast evaluation
#'   based on the available observations, \code{indRef} implements the 
#'   out-of-sample or in-sample protocol to be used and \code{generateRef} 
#'   produces the corresponding ensemble forecast given the actual observations.
#'   
#' @param nfcst number of forecast instances to be produce
#' @param type type of out-of-sample protocol to be applied (see below)
#' @param indices Subset of the observations / forecast times to be used for 
#'   reference forecasts
#' @param blocklength for cross-validation and split-sample
#'   
#' @return \item{ind}{A list of indices to be used for each forecast from 
#'   \code{1} to \code{nfcst}}
#'   
#' @section Cross-validation: Leave-one-out and leave-n-out cross-validation 
#'   reference forecasts can be produced by setting \code{type = "crossval"}. By
#'   default, the blocklength is set to \code{1}, but moving blocks of length 
#'   \code{n} can be specified by setting \code{blocklength = n}.
#'   
#' @section Split sample: In contrast to \code{type="crossval"}, 
#'   \code{type="block"} is used for split-sample validation with 
#'   non-overlapping blocks of length \code{blocklength} retained for 
#'   validation.
#'   
#' @section Forward: Correspondingly, reference forecasts that are only based on
#'   past (future) observations can be produced using \code{type = "forward"}. 
#'   For this, the first half of the reference forecasts only uses future 
#'   information, i.e. observations \code{2:n} for forecast \code{1}, \code{3:n}
#'   for \code{2} and so forth. The second half of the reference forecasts use 
#'   only past observations, i.e. observations \code{1:(n-1)} for forecast 
#'   \code{n}, \code{1:(n-2)} for \code{n-1}, etc.
#'   
#' @section Subsetting: In combination with the above, a subset of the 
#'   observations can be specified for use as reference forecasts by providing
#'   the explicit indices of the observations to be used via \code{indices=1:k}.
#'   In combination with the \code{forward} method, all observations in 
#'   \code{indices} will be used to construct the reference forecast for
#'   forecasts not included in \code{indices} (i.e. if \code{nfcst >
#'   max(indices)}).
#'   
#' @keywords utilities
#' @export
#' 
indRef <- function(nfcst, type=c('none', 'forward', 'crossval', 'block'), 
                   indices=1:nfcst, blocklength=1){
  
  ## check type of out-of-sample climatological reference generation
  type <- match.arg(type)
  stopifnot(nfcst %% 1 == 0)
  
  if (type == 'none'){
    ind <- lapply(1:nfcst, function(x) indices)
  } else if (type == 'forward') {
    stopifnot(length(indices) > 1)
    ind <- lapply(1:nfcst, function(x) indices) 
    iinds <- indices[indices %in% 1:nfcst]
    ind[iinds] <- lapply(seq(along=iinds), function(x) {
      indices[seq(ifelse(x > (length(indices) %/% 2), 1, x+1), 
                  ifelse( x > (length(indices) %/% 2), x-1, length(indices)))]
    })
  } else {
    stopifnot(blocklength < length(indices))
    if (type == 'crossval') {
      ind <- lapply(1:nfcst, function(x) setdiff(indices, 
                                                 seq(x - blocklength %/% 2,
                                                     x + (blocklength - 1) %/% 2))
      )
    } else if (type == 'block'){
      ind <- lapply(1:nfcst, function(x) indices)
      ## figure out number of blocks
      nblocks <- ceiling(nfcst / blocklength)
      for (i in seq(1, nblocks)){
        ii <- seq((i - 1)*blocklength + 1, min(nfcst, i*blocklength))
        ind[ii] <- lapply(ind[ii], function(x) setdiff(x, ii))
      }
    }
  }
  
  return(ind)
}

#' NULL
#' 
#' @param obs vector of observations
#' @param ind list or matrix of dimension (\code{n x nref}) of indices 
#'   of the observations to be used for each forecast instance
#' 
#' @export
generateRef <- function(obs, ind) {
  nobs <- length(obs)
  stopifnot(range(ind)[1] >= 1 & range(ind)[2] <= nobs)
  stopifnot(is.list(ind) | is.matrix(ind))
  
  if (is.list(ind)){
    nmax <- max(sapply(ind, length), 2)
    ind <- t(sapply(ind, function(x) c(x, rep(NA, nmax - length(x)))))
    ind <- ind[,!apply(is.na(ind), 2, all), drop=F]
  }
  
  return(array(obs[ind], dim(ind)))
}