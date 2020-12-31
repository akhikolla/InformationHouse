#' Computes Laplacian
#' 
#' @param M1 matrix of size \eqn{p \times p}
#' @param type type of laplacian (normalized (default)/combinatorial)
#' @return L Laplacian
mat.to.laplacian <- function(M1, type = "normalized"){
  
  degree.vec <- apply(M1,2,sum)
  
  if(type=="normalized"){
    du.dv.mat <- sqrt(as.matrix(degree.vec) %*% t(as.matrix(degree.vec)))
    du.dv.mat[which(du.dv.mat == 0)] <- 1
    L <- -1 * M1 / du.dv.mat
    L <- L + diag(sign(degree.vec))
  } 
  if(type=="combinatorial"){
    L <- - M1 + diag(degree.vec)
  }
  return(L)
}