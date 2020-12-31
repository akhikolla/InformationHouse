large_tensor <- function(X, Y, edges) {
  # X, Y three-dim arrays, edges a list of edges in the network
  # output is a vector res
  # in position e, corresponding to edge (z,p) via the edges-list, it is:
  # the sum over j, k of the product X_pjk * Y_zjk
  
  J <- dim(X)[2]
  K <- dim(X)[3]
  E <- nrow(edges)
  
  res <- rep(0,E)
  for (j in 1:J) {
    for (k in 1:K) {
      res <- res + X[edges[,3], j, k] * Y[edges[,2], j, k]
    }
  }
  res
}