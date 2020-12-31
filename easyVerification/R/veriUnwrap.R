# veriUnwrap.R Unwrap Arguments to Hand Over to Verification Functions
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

#' Unwrap Arguments and Hand Over to Verification Function
#' 
#' Decomposes input arguments into forecast, verifying observations, and 
#' reference forecast and hands these over to the function provided.
#' 
#' @param x n x k + 1 matrix with n forecasts of k ensemble members plus the 
#'   verifying observations
#' @param verifun character string with function name to be executed
#' @param nind named vector with number of ensemble members, ensemble members of
#'   reference forecasts, observations (defaults to 1), probability or absolute 
#'   thresholds (see details)
#' @param ref.ind list with specifications for the reference forecast
#'   (see details)
#' @param ... additional arguments passed on to \code{verifun}
#'   
#' @details Forecast verification metrics are only computed for forecasts with 
#'   non-missing verifying observation and at least one non-missing ensemble 
#'   member. Metrics for all other forecasts are set to missing. For aggregate 
#'   metrics (e.g. skill scores) the metric is computed over non-missing 
#'   observation/forecast pairs only.
#'   
#'   For computation of skill scores, reference forecasts can be provided. That 
#'   is, the first \code{nens} columns of \code{x} contain the forecasts, the 
#'   \code{(nens + 1):(ncol(x) - 1)} following columns contain the reference 
#'   forecast, and the final column contains the observations. If no reference 
#'   forecast is provided (i.e. \code{ncol(x) == nens + 1}), a climatological 
#'   forecast is constructed from the \code{n} verifying observations.
#'   
#'   The elements of vector \code{nind} have to be named with \code{nens} 
#'   containing the number of ensemble members, \code{nref} the number of 
#'   ensemble members in the reference forecast for skill scores, \code{nobs} 
#'   the number of observations (only one supported), \code{nprob} the number of
#'   probability thresholds, and \code{nthresh} the number of absolute threshold
#'   for conversion of continuous forecasts to category forecasts.
#'   
#'   \code{ref.ind} specifies the set-up of the climatological reference 
#'   forecast for skill scores if no explicit reference forecast is provided
#'   (see \code{\link{indRef}}). Also, \code{ref.ind} is used to determine the 
#'   baseline to estimate the percentile-based category boundaries to convert
#'   continuous forecasts to category probabilities.
#'      
#'   
#' @note Out-of-sample reference forecasts are now fully supported.
#'   
#' @seealso \code{\link{veriApply}}
#'   
veriUnwrap <- function(x, verifun, nind=c(nens=ncol(x) - 1, nref=0, nobs=1, nprob=0, nthresh=0), ref.ind=NULL, ...){
  nens <- nind['nens']
  nref <- nind['nref']
  nobs <- nind['nobs']
  nprob <- nind['nprob']
  nthresh <- nind['nthresh']
  stopifnot(ncol(x) >= nens + 1)
  vfun <- changearg(match.fun(verifun), format='member')
  ## extract prob or thresh
  if (nprob > 0){
    prob <- (x[,seq(nens + nref + nobs + 1,ncol(x))])[1:nprob]
  } else {
    prob <- NULL
  }
  if (nthresh > 0){
    threshold <- (x[,seq(nens + nref + nobs + 1, ncol(x))])[nprob + 1:nthresh]
  } else {
    threshold <- NULL
  }
  ## reset x
  x <- x[,seq(1,nens + nref + nobs), drop=FALSE]
  nn <- ncol(x)  
  ## mask missing values
  xmask <- apply(!is.na(x), 1, all)
  ## mask out missing values from reference construction the reference accordingly
  ref.ind <- lapply(ref.ind, function(x) cumsum(xmask)[intersect(x, which(xmask))])
  ## x <- x[xmask,,drop=FALSE] (won't work with reference forecast generation)
  ## check whether this is a skill score or a score
  is.skill <- (tolower(substr(verifun, nchar(verifun) - 1, nchar(verifun))) == 'ss' & 
                 verifun != "EnsRocss" ) |     
    substr(verifun, nchar(verifun) - 3, nchar(verifun)) == 'Diff'
  is.dress <- tolower(substr(verifun, 1, 5)) == 'dress'
  if (is.skill){
    if (nn > nens + 1){
      xref <- x[,-c(1:nens, nn), drop=F]
    } else {
      ## build reference forecast according to specifications in ref.ind
      if (is.null(ref.ind)){
        xref <- t(array(x[,nn], c(nrow(x), nrow(x))))              
      } else {
        xref <- generateRef(x[xmask,nn], ref.ind)
      }
    }
    if (is.dress){
      out <- vfun(SpecsVerification::DressEnsemble(x[xmask,1:nens,drop=FALSE]),
                  SpecsVerification::DressEnsemble(xref[xmask,,drop=FALSE]),
                  x[xmask,nn])
    } else {
      out <- vfun(convert2prob(x[xmask,1:nens,drop=FALSE], prob=prob, threshold=threshold, ref.ind=ref.ind),
                  convert2prob(xref[xmask,,drop=FALSE], prob=prob, threshold=threshold, ref.ind=ref.ind),
                  convert2prob(x[xmask,nn], prob=prob, threshold=threshold, ref.ind=ref.ind), ...)      
    }
  } else {
    stopifnot(nn == nens + 1)
    if (is.dress){
      out <- vfun(SpecsVerification::DressEnsemble(x[xmask,1:nens,drop=FALSE]),
                  x[xmask,nn])
    } else {
      out <- vfun(convert2prob(x[xmask,1:nens,drop=FALSE], prob=prob, threshold=threshold, ref.ind=ref.ind),
                  convert2prob(x[xmask,nn], prob=prob, threshold=threshold, ref.ind=ref.ind), ...)          
    }
  }
  
  ## hack to get veriUnwrap to accept non-list named vector output
  if (length(names(out)) == length(out)){
    out <- as.list(out)
  }
  

  ## check whether output has to be expanded with NA
  is.expand <- !all(xmask) & (length(out) == sum(xmask) | is.list(out))
  if (is.list(out)) is.expand <- any(sapply(out, length) == sum(xmask))
  if (is.expand){
    maskexpand <- rep(NA, length(xmask))
    maskexpand[xmask] <- 1:sum(xmask)
    if (is.list(out)){
      out <- lapply(out, function(x){
        if (length(x) == sum(xmask)){
          x <- x[maskexpand]
        }
        return(x)
      })  
    } else if (length(out) == sum(xmask)) {
      out <- out[maskexpand]
    }    
  }
  
  return(out)
}
