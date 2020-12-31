#' Construct a climatological ensemble from a vector of observations.
#'
#' @description Construct a climatological ensemble from a vector of observations. Optionally, the climatological ensemble observation at time t can be constructed without the observation at time t (leave-one-out).
#' @usage ClimEns(obs, leave.one.out=FALSE)
#' @param obs vector of length N. The observations.
#' @param leave.one.out logical, default=FALSE. If TRUE, the n-th observation is removed from the n-th row of the ensemble matrix.
#' @return matrix with N rows and N-1 columns (if leave.one.out==TRUE) or N columns otherwise. 
#' @examples
#' data(eurotempforecast)
#' ClimEns(obs)
#' @export
ClimEns <- function(obs, leave.one.out=FALSE) {

  if (length(obs) < 2 & leave.one.out == TRUE) {
    stop("Need at least 2 observations to construct leave-one-out ensemble")
  }

  # construct climatological ensemble matrix
  N <- length(obs)
  ens.clim <- t(matrix(rep(obs, N), N, N))
  
  # remove diagonal if desired
  if (leave.one.out) {
    ens.clim <- t(matrix(t(ens.clim)[-seq(1,N^2,N+1)], N-1, N))
  }
  
  return(ens.clim)
}

