#' Updates connection sign matrix.
#' @param M penalty matrix
#' @param xi estimated \eqn{\hat{\xi}_{ij}} matrix
#' @return M updated correct M matrix
matrix.M.update <- function(M,xi){
  init       <- -(xi)*abs(M) # update M1 matrix
  diag(init) <- abs(diag(M)) # store diagonal values as in initial penalty matrix
  return(init)
}