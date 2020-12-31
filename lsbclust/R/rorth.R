#' Generate A Random Orthonormal Matrix
#' 
#' Uniformly sample an orthornormal matrix from the collection of all possible 
#' orthonormal matrices of a certain size. The QR decomposition is used on a 
#' matrix containing Gaussian random numbers. The QR decomposition might not
#' be the most efficient algorithm under some circumstances.
#' 
#' @references 
#' Stewart, G. W. (1980). The efficient generation of random orthogonal matrices with an 
#' application to condition estimators. \emph{SIAM Journal on Numerical Analysis}, 17(3), 403-409.
#' 
#' @param nrow Integer giving the number of rows required.
#' @param ncol Integer giving the number of columns required.
#' @param sd The standard deviation passed to \code{\link{rnorm}}
#' @export
#' 
#' @examples
#' set.seed(1)
#' rorth(5, 2)
rorth <- function(nrow, ncol, sd = 1) {
  
  ## Generate matrix, compute QR decomposition
  mat <- matrix(rnorm(nrow * ncol, sd = sd), nrow = nrow, ncol = ncol)
  mat <- qr(mat)
  Q <- qr.Q(mat)
  R <- qr.R(mat)
  
  ## Assure that diagonal of R is positive
  Q <- Q %*% diag(sign(diag(R)))
  
  return(Q)
}