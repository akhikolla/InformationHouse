##' @title Empirical Highest Posterior Density Interval
##' @description Computes the Highest Posterior Density (HPD) interval from a posterior sample. Works only for scalar marginal posteriors.
##' @usage hpd(x, prob = 0.95)
##'
##' @param x a univariate or a posterior sample for a scalar parameter.
##' @param prob posterior probability content.
##'
##' @return vector of two doubles.
##' @export
hpd <- function(x, prob = 0.95){
  conf <- min(prob, 1-prob)
  n <- length(x)
  nn <- round( n*conf )
  x <- sort(x)
  xx <- x[ (n-nn+1):n ] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx==m)[1]
  return( c( x[ nnn ], x[ n-nn+nnn ] ) )
}
