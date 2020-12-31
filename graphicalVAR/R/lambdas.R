invGlasso <- function(x){
  if (all(eigen(x)$values > sqrt(.Machine$double.eps))){
    Xinv <- solve(x)
  } else {
    Xglas <- glasso(x,0.05,penalize.diagonal=FALSE)    
    Xinv <- Xglas$wi
  }
  Xinv
}
generate_lambdas <- function(
  X,
  Y,
  nLambda_kappa = 10,
  nLambda_beta = 10,
  lambda_min_kappa = 0.05,
  lambda_min_beta = 0.05, 
  penalize.diagonal=TRUE,
  version0.1.4 = FALSE
){
  N <- nrow(Y)
  P <- ncol(Y)
  

  #### Lambda sequence for Kappa:
  corY <- cov2cor(t(Y)%*%Y/nrow(Y))
  if (version0.1.4){
    lam_K_max = max(abs(corY))
  } else {
    lam_K_max = max(abs(corY[upper.tri(corY)]))
  }
  lam_K_min = lambda_min_kappa*lam_K_max
  lam_K = exp(seq(log(lam_K_max), log(lam_K_min), length = nLambda_kappa))
  
  #### Lambda sequence for Beta
  # Initial estimate for Kappa:
  # Yinv <- invGlasso(t(Y) %*% Y / N)
  # Xinv <- invGlasso(t(X) %*% X / N)
#   beta <- t(Y) %*% X %*% Xinv
#   S <- 1/(nrow(Y)) * (
#     t(Y) %*% Y -
#       t(Y) %*% X %*% t(beta) -
#       beta %*% t(X) %*% Y +
#       beta %*% t(X) %*% X %*% t(beta)
#   )
#   S <- (S + t(S)) / 2
#   if (any(eigen(S)$value < -sqrt(.Machine$double.eps))) stop("Residual covariances not postive definite")
#   kappa <- invGlasso(S)
#   kappa <- (kappa + t(kappa)) / 2
  # lam_B_max = max(abs((1/N)*t(X)%*%Y%*%Yinv))
  
  Yinv <- invGlasso(t(Y) %*% Y)
  lam_B_max = max(abs(t(X)%*%Y%*%Yinv))
  lam_B_min = lambda_min_beta*lam_B_max
  lam_B = exp(seq(log(lam_B_max), log(lam_B_min), length = nLambda_beta))
  
  return(list(lambda_kappa = lam_K, lambda_beta = lam_B))
}
  

### BACKWARD COMPETABILITY ###

SparseTSCGM_lambdas <-
  function(X, Y, nlambda = 100, lambda.min.ratio = 0.01){
    # From SparseTSCGM package:
    lambda.seq <- function(SS, SA,nlambda)
    {
      if (length(nlambda)==1) nlambda <- rep(nlambda,2)
      # lambda.min.ratio=0.1
      d= dim(SS)[2]
      lambda.max1 = max(max(SS-diag(d)),-min(SS-diag(d)))
      lambda.min1 = lambda.min.ratio*lambda.max1
      lambda1 = exp(seq(log(lambda.max1), log(lambda.min1), length = nlambda[1]))
      lambda.min.ratio2=0.15
      lambda.max2 = max(max(SA),-min(SA))
      lambda.min2 = lambda.min.ratio2*lambda.max2
      lambda2 = exp(seq(log(lambda.max2), log(lambda.min2), length = nlambda[2]))
      return(list(lambda1=lambda1, lambda2=lambda2))
    }
    
    T <- dim(Y)[1]
    p <- dim(X)[2]
    n <- 1
    q <- dim(Y)[2]
    xtyi <- array(NA, c(p,q,n))
    xtxi <- array(NA, c(p,p,n))
    ytyi <- array(NA, c(q,q,n))
    
    XX <- X
    YY <- Y
    XX2 <- X^2
    YY2 <- Y^2
    xtyi <- crossprod(XX,YY)
    xtxi <- crossprod(XX)
    ytyi <- crossprod(YY)
    
    xty=apply(xtyi, c(1,2), sum)
    xtx=apply(xtxi, c(1,2), sum)
    yty=apply(ytyi, c(1,2), sum)
    xtxt=apply(xtxi, c(1,2), sum)/(n*T)
    xtx2=(n*T)*colMeans(apply(XX2, c(1,2), sum))
    yty2=(n*T)*colMeans(apply(YY2, c(1,2), sum))
    
    SX <- xtx/(n*T)
    mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
    
    SX <- xtx/(n*T)
    mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
    SXi <- mSX$wi
    SS =(yty)/(n*T)
    SS = cov2cor(SS)
    SAs = xty/(n*T)
    
    SA = t(SAs) %*% SXi
    
    lambda <-  lambda.seq(SS=SS,SA=SA, nlambda=nlambda)
    lam1 <- round(lambda$lambda1,3) 
    lam2 <- round(lambda$lambda2,3)
    lam2 <- round(lam2/max(lam2),3)    
    return(list(lambda_kappa = lam1, lambda_beta = lam2))
  }
  
  