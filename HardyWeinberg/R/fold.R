fold <-
function(X,lower=TRUE) {
  nr <- nrow(X)
  nc <- ncol(X)
  if(nr!=nc) {
    stop("X is not square.")
  } else {
    lo <- X[lower.tri(X, diag = FALSE)]
    Xt <- t(X)
    up <- Xt[lower.tri(Xt, diag = FALSE)]
    X[lower.tri(X, diag = FALSE)] <- lo+up
    X[upper.tri(X, diag = FALSE)] <- 0
  }
  if(!lower) X <- t(X) 
  return(X)
}
