#' Randomly Generate Positive Singular Values
#' 
#' Generate random singular values for a specified number of clusters
#' for use in simulations. A mixture distribution is used with truncation to ensure that the
#' singular values differ between clusters, are ordered, and are nonnegative.
#' 
#' @param nclust Integer giving the number of clusters for which to sample singular values.
#' @param ndim Integer; the number of singular values required.
#' @param mins Numeric vector of length \code{ndim} giving the minimum values for the 
#' respective singular values.
#' @param max Numeric value giving the maximum possible value for the mean of the cluster-
#' specific singular value distribution, relative to the \code{mins}
#' @export
simsv <- function(nclust, ndim = 2, mins = 1, max = 5) {
  
  ## Expand mins if needed
  if (length(mins) == 1)
    mins <- rep(mins, ndim)
  
  ## Generate the means
  mns <- matrix(runif(nclust * ndim, max = max), nrow = nclust, ncol = ndim)
  mns <- t(apply(mns, 1, sort, decreasing = TRUE))
  mns <- rep(1L, nclust) %o% mins + mns
  # plot(mns, asp = 1)
  # abline(a = mins[2], b = 1)
  
  ## Return
  return(split(mns, seq_len(nclust)))
}