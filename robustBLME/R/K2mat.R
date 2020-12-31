
##' @title The K2 (n x n) matrix
##' @description Computes the K2 matrix that is needed for the consistency
##' correction of the variance prameters in the robust REML2 estimating equations.
##' @usage K2mat(c2_hub, n)
##' @param c2_hub Robustness constant for the Huber's proposal 2 estimating equation.
##' @param n The number of rows and columns of the matrix.
##' @return matrix
##' @noRd
K2mat <- function(c2_hub, n){
  ss <- integrate(f = Vectorize(FUN = function(x) dnorm(x)*psi_huber(x, c = c2_hub) ^ (2)),
                  lower = -10, upper = 10)$value
  st = 0.0
  ans = matrix(st, n, n)
  diag(ans) = ss
  return(ans)
}