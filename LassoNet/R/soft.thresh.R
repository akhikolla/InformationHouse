#' Evaluates each element of \eqn{\beta} using soft thresholding operator.
#' @param x \eqn{\beta} coordinate
#' @param kappa \eqn{\kappa} value in general or \eqn{\lambda_1} for covariance updating.
#' @return x value after applying soft thresholding operator
soft.thresh <- function(x ,kappa){
  # initializing
  
  # changing according to rule
  if (abs(x)>kappa){
    x <- sign(x)*(abs(x)-kappa)
  } else{
    x <- 0
  }
  
  return(x)
}