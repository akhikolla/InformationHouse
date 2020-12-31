

randomGVARmodel <- function(
  Nvar,
  probKappaEdge = 0.1,
  probKappaPositive = 0.5,
  probBetaEdge = 0.1,
  probBetaPositive = 0.5,
  maxtry = 10,
  kappaConstant = 1.1
){
  try <- 0
  
  repeat{
    kappaRange = c(0.5,1)
    ## Approach from 
    # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.
    trueKappa <- matrix(0,Nvar,Nvar)
    trueKappa[upper.tri(trueKappa)] <- sample(c(0,1),sum(upper.tri(trueKappa)),TRUE,prob=c(1-probKappaEdge,probKappaEdge))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(1,-1),sum(upper.tri(trueKappa)),TRUE,prob=c(1-probKappaPositive,probKappaPositive)) * 
      runif(sum(upper.tri(trueKappa)), min(kappaRange),max(kappaRange))
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]  
    diag(trueKappa) <- kappaConstant * rowSums(abs(trueKappa))
    diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
    trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
    trueKappa <- (trueKappa + t(trueKappa)) / 2
    
    Vmin <- min(abs(trueKappa[trueKappa!=0]))
    trueBeta <- matrix(sample(c(0,1),Nvar^2,TRUE,prob=c(1-probBetaEdge,probBetaEdge)),Nvar,Nvar)
    trueBeta <- trueBeta * sample(c(-1,1),Nvar^2,TRUE,prob=c(1-probBetaPositive,probBetaPositive)) * 
      runif(Nvar^2, Vmin,1)
    
    diag(trueBeta) <- Vmin
    
    evK <- eigen(trueKappa)$values
    evB <- eigen(trueBeta)$values
    
    while(any(Re(evB)^2 + Im(evB)^2 > 1)){
      trueBeta <- 0.95*trueBeta
      evB <- eigen(trueBeta)$values
    }
    
    
    if (all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){
      break
    }
    
    try <- try + 1
    if (try > maxtry){
      stop("Maximum number of tries reached.")
    }
    
  }
  Res <- list(
    kappa = trueKappa,
    beta = trueBeta,
    PCC = computePCC(trueKappa),
    PDC = computePDC(trueBeta,trueKappa)
  )
  class(Res) <- "gVARmodel"
  return(Res)
}



graphicalVARsim <- function(
  nTime, # Number of time points
  beta, # if dim  is 2xnVar, assume changescores
  kappa,
  mean = rep(0,ncol(kappa)),
  # sd = rep(1, ncol(kappa)),
  init = mean,
  warmup = 100,
  lbound = rep(-Inf, ncol(kappa)),
  ubound = rep(Inf, ncol(kappa))
#   skewed = rep(0,ncol(kappa)),
#   WN = FALSE
){
  
  
  stopifnot(!missing(beta))
  stopifnot(!missing(kappa))  
  
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  
  totTime <- nTime + warmup
  
  Data <- t(matrix(init, Nvar, totTime))
  
  Sigma <- solve(kappa)
  
  #   lbound <- (lbound - mean) / sd
  #   ubound <- (ubound - mean) / sd
  #   
  
  # if (sum(skewed==rep(0,ncol(kappa)))==ncol(kappa)){
    for (t in 2:totTime){
      Data[t,] <- mean + t(beta %*% (Data[t-1,]-mean))  + mvtnorm::rmvnorm(1, rep(0,Nvar), Sigma)
      Data[t,] <- ifelse(Data[t,]  < lbound, lbound, Data[t,] )
      Data[t,] <- ifelse(Data[t,]  > ubound, ubound, Data[t,] )
    }
#   }else{
#     for (t in 2:totTime){#Needed to round Omega to avoid error "not symmetrical"
#       Data[t,] <- t(beta %*% Data[t-1,])  + sn::rmsn(n=1,dim=1*Nvar, mu=rep(0,Nvar),Omega=round(Sigma,digits = 5),
#                                                      alpha=skewed)[1,]
#       Data[t,] <- ifelse(Data[t,]  < lbound, lbound, Data[t,] )
#       Data[t,] <- ifelse(Data[t,]  > ubound, ubound, Data[t,] )
#     }
#   }
#   if (WN){
#     for (t in 2:totTime){
#       Data[t,] <- Data[t,] +  rmvnorm(1, rep(0,Nvar), Sigma)
#     }
#   }
 
  
  return(Data[-seq_len(warmup), ,drop=FALSE])
}

# 
# graphicalVARsim <- function(
#   nTime, # Number of time points
#   beta, # if dim  is 2xnVar, assume changescores
#   kappa,
#   scaledMatrices = TRUE,
#   mean = rep(0,ncol(kappa)),
#   sd = rep(1, ncol(kappa)),
#   init = mean,
#   intercepts = 0,
#   warmup = 100,
#   pi =rep(0,ncol(kappa)),
#   nudgeMean = rep(0.1,ncol(kappa)),
#   nudgeSD = rep(0.1,ncol(kappa)),
#   lbound = rep(-Inf, ncol(kappa)),
#   ubound = rep(Inf, ncol(kappa))
# ){
#   
#   
#   stopifnot(!missing(beta))
#   stopifnot(!missing(kappa))  
#   
#   Nvar <- ncol(kappa)
#   init <- rep(init, length = Nvar)
#   intercepts <- rep(intercepts, length = Nvar)
#   
#   
#   # Temp solulution, add change scores:
#   if (ncol(beta)==Nvar){
#     beta <- cbind(beta,matrix(0,Nvar,Nvar))
#   }
#   if (ncol(beta)!=Nvar*2) stop("Beta must contain only lag-1 and changescore effects.")
#   
#   totTime <- nTime + warmup
#   
#   Data <- t(matrix(init, Nvar, totTime))
#   
#   Sigma <- solve(kappa)
#   
#   lbound <- (lbound - mean) / sd
#   ubound <- (ubound - mean) / sd
#   
#   for (t in 3:totTime){
#     residMean <- ifelse(runif(Nvar) < pi, rnorm(Nvar,nudgeMean,nudgeSD), 0)
#     
#     Data[t,] <- t(intercepts + beta[,1:Nvar] %*% Data[t-1,] + beta[,Nvar + (1:Nvar)] %*% (Data[t-1,]-Data[t-2,])) + rmvnorm(1, residMean, Sigma)
#     Data[t,] <- ifelse(Data[t,]  < lbound, lbound, Data[t,] )
#     Data[t,] <- ifelse(Data[t,]  > ubound, ubound, Data[t,] )
#   }
#   
#   for (t in 1:totTime){
#     Data[t,] <- Data[t,]*sd + mean
#   }
#   
#   return(Data[-seq_len(warmup), ,drop=FALSE])
# }