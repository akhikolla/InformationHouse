toTriangular <- function(x) {
  x <- unlist(x)
  k <- n.alleles(x)
  Ma <- matrix(rep(0,k*k),ncol=k)
  colnames(Ma) <- LETTERS[1:k]
  rownames(Ma) <- colnames(Ma)
  Ma[lower.tri(Ma, diag=TRUE)] <- x
  return(Ma)
}
