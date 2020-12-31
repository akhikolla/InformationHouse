#' Compute least squares coefficient fast
#' 
#' @params y response vector
#' @params x matrix of covariates
#' @return param parameter vector
fastols <- function (y ,x) { 
  XtX <- crossprod(x) 
  Xty <- crossprod(x, y) 
  return(solve(XtX, Xty)) 
}
