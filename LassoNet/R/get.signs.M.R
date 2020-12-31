#' Computes vector form of connection sign matrix for storage purposes
#' 
#' @param MAT A matrix of connection signs that contains -1, 1 or 0
#' @return vec.out Vectorized MAT matrix
get.signs.M <- function(MAT) {
  vec.out <- cbind(which(MAT!=0, arr.ind = T),as.matrix(MAT[which(MAT!=0, arr.ind = T)]))
  return(vec.out)
}
