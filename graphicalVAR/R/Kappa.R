Kappa <-
function(beta, X, Y, lambda_kappa){
  n <- nrow(Y)  
  SigmaR <- 1/n * t(Y - X %*% beta) %*% (Y - X %*% beta)
  if (any(eigen(SigmaR,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    stop("Residual covariance matrix is not non-negative definite")
  }
  res <- glasso(SigmaR, lambda_kappa, penalize.diagonal = FALSE)
  return(as.matrix(forceSymmetric(res$wi)))
}
