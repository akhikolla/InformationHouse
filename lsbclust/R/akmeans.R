#' K-Means Over One Way of An Three-Way Array
#' 
#' Vectorize matrix slices over a specific way of an three-way array, and
#' conduct \code{\link{kmeans}} on it.
#' 
#' 
#' @param data Three-way data array
#' @param margin Integer indicating which way to cluster over
#' @param centers Passed to \code{\link{kmeans}}
#' @param ndim The rank of the low dimensional approximation of the matrix
#' slices to construct before clustering (using \code{\link{svd}})
#' @param \dots Additional arguments passed to \code{\link{kmeans}}
#' 
#' @export
#' @importFrom plyr aaply
#' @examples 
#' set.seed(1)
#' res <- akmeans(data = carray(dcars), margin = 3L, centers = 5, nstart = 10)
akmeans <- function(data, centers, margin = 3L, ndim = NULL, ...) {
  
  ## Checks
  if (!is(data, "array") || length(dim(data)) != 3L)
    stop("Argument 'data' must be a three-way array.")
  
  ## Preliminaries
  dims <- dim(data)
  
  ## Do low-rank decomposition
  if (!is.null(ndim)) {
    approxfun <- function(x) {
      sv <- svd(x)
      return(sv$u[, 1:ndim, drop = FALSE] %*% diag(sv$d[1:ndim], nrow = ndim, ncol = ndim) %*% 
               t(sv$v[, 1:ndim, drop = FALSE]))
    }
    data <- aaply(.data = data, .margins = margin, .fun = approxfun)
    data <- aperm(data, perm = order(c(margin, (1:3)[-margin])))
  }
  
  ## Convert data to matrix by vectorising the slices
  # x <- matrix(data, nrow = dims[margin], ncol = prod(dims[-margin]), byrow = TRUE)
  x <- t(apply(data, MARGIN = margin, FUN = identity))
  res <- kmeans(x = x, centers = centers, ...)
  
  ## Get fitted and convert back to array
  fitted <- fitted(res)
  fitted <- array(fitted, dim = c(dims[margin], dims[-margin]))
  fitted <- aperm(fitted, perm = order(c(margin, (1:3)[-margin])))
  res$fitted <- fitted
  
  ## Convert centers back to an array
  centers <- array(t(res$centers), dim = c(dims[-margin], centers))
  
  ## Permute such that margin is correct
  centers <- aperm(centers, perm = order(c((1:3)[-margin], margin)))
  res$centers_array <- centers
  
  
  ## Return
  class(res) <- c("akmeans", class(res))
  return(res)
}