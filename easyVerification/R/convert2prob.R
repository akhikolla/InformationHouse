# convert2prob.R Convert to Category Forecast
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

#' @name convert2prob
#' 
#' @aliases prob2thresh
#' 
#' @title Convert to Probability / Category Forecast
#' 
#' @description \code{convert2prob} Converts the continuous ensemble forecast to
#' counts of ensemble members per category. The categories can be defined
#' relative to the ensemble distribution (using \code{prob}) or relative to
#' absolute values for the category thresholds (using \code{threshold}, see
#' details). \code{prob2thresh} converts the relative threshold to absolute
#' thresholds for later processing. \code{expandthresh} expands the vector or
#' matrix of thresholds to fit the input data.
#' 
#' @param x input vector or matrix
#' @param prob thresholds for categorical forecasts (defaults to NULL)
#' @param threshold absolute thresholds for categorical forecasts (defaults to 
#'   NULL)
#' @param ref.ind list of forecast/obs instances to be used to estimate 
#'   percentile thresholds
#' @param multi.model logical, are we dealing with initial condition (the 
#'   default) or multi-model ensembles (see details)?
#'   
#' @details In case both \code{prob} and \code{threshold} are set to 
#'   \code{NULL}, the function returns the input \code{x} without modification.
#'   If \code{prob} is set, a matrix with the number of occurrences per class for
#'   a given quantile of the full distribution (e.g. temperature above/below the
#'   median). If \code{threshold} is set, the classes are defined based on the 
#'   absolute value (e.g. temperature above/below 13 deg. C). Multiple classes
#'   are
#'   
#'   Only certain formats of \code{threshold} and \code{prob} are supported.
#'   \code{prob} has to be a vector with percentile thresholds separating the
#'   different classes. \code{threshold} can be a vector, matrix or array with
#'   the first entry corresponding to the different classes, and the last to the
#'   different ensemble members (if present). Thereby, time/forecast varying
#'   thresholds can potentially be supplied (although I am not sure this is
#'   useful or needed).
#'   
#'   If \code{ref.ind} is specified, only the specified indices of the input
#'   variables are used to estimate the percentile thresholds (\code{prob}). If
#'   used with \code{threshold}, or without anything, \code{ref.ind} has no effect. 
#'   
#'   If \code{multi.model = TRUE}, the relative thresholds supplied by
#'   \code{prob} are ensemble member specific, i.e. are estimated for each
#'   ensemble member separately. This is in particular applicable for
#'   multi-model ensembles with model dependent biases.
#'   
#' @return Matrix of occurrences per class (i.e. the number of ensemble members 
#'   per class, or an indicator for the observations)
#'   
#' @examples
#' tm <- toymodel()
#' 
#' ## convert to tercile forecasts (only display first forecast and obs)
#' convert2prob(tm$fcst, prob=1:2/3)[1,]
#' convert2prob(tm$obs, prob=1:2/3)[1,]
#' 
#' ## convert to category forecasts (smaller and larger than 1)
#' convert2prob(tm$fcst, threshold=1)[1,]
#' convert2prob(tm$obs, threshold=1)[1,]
#' 
#' @seealso \code{\link{veriApply}}
#'   
#' @keywords utilities
#' @export
convert2prob <- function(x, prob=NULL, threshold=NULL, ref.ind=NULL,
                         multi.model=FALSE){
  ## check if input arguments comply
  stopifnot(is.vector(x) | is.matrix(x))
  stopifnot(any(!is.na(x)))

  ## do nothing if you do not need to do something, simple things first
  if (is.null(prob) & is.null(threshold)) return(x)
  if (!is.null(prob) & !is.null(threshold)){
    stop('Both probability and absolute thresholds provided')
  } 
  
  if (!is.null(prob)){
    ## convert probability to absolute threshold
    stopifnot(unlist(ref.ind) %in% 1:nrow(as.matrix(x)))
    threshold <- prob2thresh(x=x, prob=prob, ref.ind=ref.ind, multi.model=multi.model)
  } else {
    ## blow up  threshold to size of nclass x size(x)
    if (is.null(prob)) threshold <- expandthresh(threshold, x)
  }

  ## count occurrences in bins
  nclass <- nrow(threshold) + 1
  xtmp <- array(apply(rep(x, each=nrow(threshold)) > threshold, -1, sum), dim(as.matrix(x))) + 1
  xout <- t(apply(xtmp, 1, tabulate, nbins=nclass))      
  xout[apply(as.matrix(is.na(x)), 1, any),] <- NA
  
  return(xout)
}


#' @rdname convert2prob
prob2thresh <- function(x, prob, ref.ind=NULL, multi.model=FALSE){
  
  ## reduce size of x if constructed from observations (i.e. all members the same)
  ## to guarantee a consistent estimate of the quantiles
  xthresh <- x
  if (is.matrix(x)){
    if (is.null(ref.ind)){
      if(all(apply(x, 1, function(y) all(y == x[1,])))){
        xthresh <- x[1,]
      }
    } else {
      maxind <- max(unlist(ref.ind))
      if (length(unique(c(x))) <= maxind) {
        ## generate indices of reference forecast to reverse
        ## engineer original vector
        iref <- generateRef(seq(1, maxind), ref.ind)
        if (all(dim(iref) == dim(x))){
          if (all(tapply(x, iref, function(y) all(y == y[1])))){
            iind <- which(!duplicated(c(iref)))
            xthresh <- x[iind][order(iref[iind])]
          }
        }
      }
    }
  }
  
  ## apply out-of-sample strategy if needed
  ## already adjust the size of threshold to size of output to minimize errors due to ambiguity
  if (length(unique(ref.ind)) <= 1){ ## condition for both is.null or all equal
    if (is.null(ref.ind)) ref.ind <- list(1:nrow(as.matrix(xthresh)))
    if (multi.model){
      threshold <- apply(as.matrix(xthresh)[ref.ind[[1]],,drop=F], 2, quantile, sort(prob), na.rm=T, type=8)
      threshold <- array(threshold[rep(1:nrow(threshold), length=length(prob)*nrow(as.matrix(x))),], 
                         c(length(prob), size(x)))
    } else {
      threshold <- array(quantile(as.matrix(xthresh)[ref.ind[[1]],], sort(prob), na.rm=T, type=8), 
                         c(length(prob), size(x))) 
    }
  } else {
    if (multi.model){
      threshold <- aperm(sapply(seq(along=ref.ind), function(i){
        apply(as.matrix(xthresh)[ref.ind[[i]],, drop=FALSE], 2, quantile, sort(prob), na.rm=T, type=8)
      }, simplify='array'), c(1,3,2))
    } else {
      threshold <- array(sapply(seq(along=ref.ind), function(i){
        quantile(as.matrix(xthresh)[ref.ind[[i]],], sort(prob), na.rm=T, type=8)
      }, simplify='array'), c(length(prob), size(x)))
    }
  }
  return(threshold)
}

#' @rdname convert2prob
expandthresh <- function(threshold, x){
  nclass <- size(threshold)[1]
  if (is.vector(threshold)){
    threshold <- array(threshold, c(nclass, size(x)))
  } else if (! all(size(threshold)[-1] == size(x))){
    stopifnot(ncol(threshold) == ncol(as.matrix(x)))
    threshold <- array(threshold[,rep(1:ncol(threshold), each=nrow(x))],
                       c(nrow(threshold), size(x)))
  } 
  return(threshold)
}
