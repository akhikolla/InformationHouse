Rothmana <-
function(X, Y, lambda_beta, lambda_kappa, convergence = 1e-4, gamma = 0.5, maxit.in = 100, maxit.out = 100,
         penalize.diagonal, # if FALSE, penalizes the first diagonal (assumed to be auto regressions), even when ncol(X) != ncol(Y) !
         interceptColumn = 1, # Set to NULL or NA to omit
         mimic = "current",
         likelihood = c("unpenalized","penalized")
         ){
  # Algorithm 2 of Rothmana, Levinaa & Ji Zhua
  
  likelihood <- match.arg(likelihood)
  
  nY <- ncol(Y)
  nX <- ncol(X)
 
  if (missing(penalize.diagonal)){
    if (mimic == "0.1.2"){
      penalize.diagonal <- nY != nX
    } else {
      penalize.diagonal <- (nY != nX-1) & (nY != nX ) 
    }
  }
  
  lambda_mat <- matrix(lambda_beta,nX, nY)
  if (!penalize.diagonal){
    if (nY == nX){
      add <- 0
    } else if (nY == nX - 1){
      add <- 1
    } else {
      stop("Beta is not P x P or P x P+1, cannot detect diagonal.")
    }
    for (i in 1:min(c(nY,nX))){
      lambda_mat[i+add,i] <- 0
    }
  }
  if (!is.null(interceptColumn) && !is.na(interceptColumn)){
    lambda_mat[interceptColumn,] <- 0
  }
 
  n <- nrow(X)
  beta_ridge <- beta_ridge_C(X, Y, lambda_beta)
  
  # Starting values:
  beta <- matrix(0, nX, nY)  
  
  # Algorithm:
  it <- 0

  repeat{
    it <- it + 1
    kappa <- Kappa(beta, X, Y, lambda_kappa)
    beta_old <- beta
    beta <- Beta_C(kappa, beta, X, Y, lambda_beta, lambda_mat, convergence, maxit.in) 
    
    if (sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){
      break
    }
    
    if (it > maxit.out){
      warning("Model did NOT converge in outer loop")
      break
    }
  }
  
  ## Compute unconstrained kappa (codes from SparseTSCGM):
  ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
  WS <-  (t(Y)%*%Y - t(Y) %*% X  %*% beta - t(beta) %*% t(X)%*%Y + t(beta) %*% t(X)%*%X %*% beta)/(nrow(X))
  
  if (any(eigen(WS,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    stop("Residual covariance matrix is not non-negative definite")
  }
  
  if (likelihood == "unpenalized"){
    if (nrow(ZeroIndex)==0){
      out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
    } else {
      out4 <- suppressWarnings(glasso(WS, rho = 0, zero = ZeroIndex,
                                      trace = FALSE))
    }
    lik1  <- determinant( out4$wi)$modulus[1]
    lik2 <- sum(diag( out4$wi%*%WS))
  } else {
    lik1  <- determinant( kappa )$modulus[1]
    lik2 <- sum(diag( kappa%*%WS))
  }

  pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
  if (mimic == "0.1.2"){
    pdB = sum(sum(beta !=0))
  } else {
    pdB = sum(sum(beta[lambda_mat!=0] !=0)) 
  }
  
  LLk <-  (n/2)*(lik1-lik2) 
  LLk0 <-  (n/2)*(-lik2)
  
  EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*nY)

  
  ### TRANSPOSE BETA!!!
  return(list(beta=t(beta), kappa=kappa, EBIC = EBIC))
}
