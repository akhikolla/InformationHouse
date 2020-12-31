#' Calculate predictive weights for CGGP
#' 
#' Predictive weights are Sigma^{-1}*y in standard GP.
#' This calculation is much faster since we don't need to
#' solve the full system of equations.
#'
#' @param CGGP CGGP object
#' @param y Measured values for CGGP$design
#' @param theta Correlation parameters
#' @param return_lS Should lS be returned?
#'
#' @return Vector with predictive weights
#' @export
#'
#' @examples
#' cggp <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(cggp$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' CGGP_internal_calcpw(CGGP=cggp, y=y, theta=cggp$thetaMAP)
CGGP_internal_calcpw <- function(CGGP, y, theta, return_lS=FALSE) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max value of all blocks
  # Now going to store choleskys instead of inverses for stability
  #CiS = list(matrix(1,1,1),Q*CGGP$d) # A list of matrices, Q for each dimension
  
  cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
  # Loop over each dimension
  for (dimlcv in 1:CGGP$d) {
    # Loop over each possible needed correlation matrix
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = FALSE)
      S = Sstuff
      # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
      solvetry <- try({
        cS = chol(S)
        cholS[[(dimlcv-1)*Q+levellcv]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      }, silent = TRUE)
      if (inherits(solvetry, "try-error")) {return(Inf)}
      lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
    }
  }
  
  if(!is.matrix(y)){
    pw = rep(0, length(y)) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        B = y[IS]
        rcpp_kronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
        pw[IS] = pw[IS]+CGGP$w[blocklcv] * B
      }
    }
    if (return_lS) {
      return(list(pw=pw, lS=lS))
    }else{
      return(pw)
    }
  }else{
    numout = dim(y)[2]
    pw = matrix(0,nrow=dim(y)[1],ncol=numout) # Predictive weight for each measured point
    # Loop over blocks selected
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        VVV1 = unlist(cholS[gg+CGGP$uo[blocklcv,]]);
        VVV2 = CGGP$gridsizest[blocklcv,];
        for(outdimlcv in 1:numout){
          B = y[IS,outdimlcv]
          rcpp_kronDBS(VVV1, B, VVV2)
          pw[IS,outdimlcv] = pw[IS,outdimlcv]+CGGP$w[blocklcv] * B
        }
      }
    }
    if (return_lS) {
      return(list(pw=pw, lS=lS))
    }else{
      return(pw)
    }
  }
}

#' Calculate derivative of pw
#'
#' @inheritParams CGGP_internal_calcpw
#' @param return_lS Should lS and dlS be returned?
#'
#' @return derivative matrix of pw with respect to logtheta
#' @export
#' @import Rcpp
#'
#' @examples
#' cggp <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(cggp$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' CGGP_internal_calcpwanddpw(CGGP=cggp, y=y, theta=cggp$thetaMAP)
CGGP_internal_calcpwanddpw <- function(CGGP, y, theta, return_lS=FALSE) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*CGGP$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$numpara*CGGP$d)
  }
  
  # Loop over each dimension
  for (dimlcv in 1:CGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:CGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        for(paralcv in 1:CGGP$numpara){
          dlS[levellcv, CGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
        }
      }
    }
  }
  
  pw = rep(0, length(y)) # predictive weights
  dpw = matrix(0, nrow = CGGP$numpara*CGGP$d, ncol = length(y)) # derivative of predictive weights
  gg = (1:CGGP$d-1)*Q
  for (blocklcv in 1:CGGP$uoCOUNT) {
    if(abs(CGGP$w[blocklcv])>0.5){
      IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
      B = CGGP$w[blocklcv]*y[IS]
      dB = rcpp_gkronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
      dpw[,IS] = dpw[,IS] +dB
      pw[IS] = pw[IS] + B
    }
  }
  dpw =t(dpw)
  out <- list(pw=pw,
              dpw=dpw)
  if (return_lS) {
    out$lS <- lS
    out$dlS <- dlS
  }
  
  out
}



CGGP_internal_calcsigma2 <- function(CGGP, y, theta, return_lS=FALSE) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  if(return_lS){
    lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
  }
  
  # Loop over each dimension
  for (dimlcv in 1:CGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = FALSE)
      S = Sstuff
      # cS = chol(S)
      cS = try(chol(S))
      if (inherits(cS, "try-error")) {
        stop("Cholesky error in CGGP_internal_calcsigma2")
      }
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
      }
    }
  }
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+CGGP$uo[blocklcv,]])
        VVV3=CGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          B0 = y[IS,outdimlcv]
          B = (CGGP$w[blocklcv]/dim(y)[1])*B0
          rcpp_kronDBS(VVV1,B,VVV3)
          sigma2[outdimlcv] = sigma2[outdimlcv]  + t(B0)%*%B
        }
      }
    }
    out <- list(sigma2=sigma2)
    if (return_lS) {
      out$lS <- lS
    }
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=CGGP$d) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        B0 = y[IS]
        B = (CGGP$w[blocklcv]/length(y))*B0
        rcpp_kronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]),B, CGGP$gridsizest[blocklcv,])
        sigma2 = sigma2 + t(B0)%*%B
        if (any(is.na(sigma2))) {warning("sigma2 is NA in CGGP_internal_calcsigma2")}
      }
    }
    out <- list(sigma2=sigma2)
    if (return_lS) {
      out$lS <- lS
    }
  }
  return(out)
}

CGGP_internal_calcsigma2anddsigma2 <- function(CGGP, y, theta, return_lS=FALSE) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max level of all blocks
  cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
  dMatdtheta = list(matrix(1,1,1),Q*CGGP$d)
  
  if(return_lS){
    lS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$d) # Save log determinant of matrices
    dlS = matrix(0, nrow = max(CGGP$uo[1:CGGP$uoCOUNT,]), ncol = CGGP$numpara*CGGP$d) 
  }
  
  # Loop over each dimension
  for (dimlcv in 1:CGGP$d) {
    # Loop over depth of each dim
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
      Xbrn = Xbrn[order(Xbrn)]
      nv = length(Xbrn);
      Sstuff = CGGP$CorrMat(Xbrn, Xbrn , theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = TRUE)
      S = Sstuff$C
      
      cS = chol(S)
      
      cholS[[(dimlcv-1)*Q+levellcv]] = cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
      dMatdtheta[[(dimlcv-1)*Q+levellcv]] = -backsolve(cS,backsolve(cS,Sstuff$dCdtheta, transpose = TRUE))
      for(paralcv in 1:CGGP$numpara){
        dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv] = t(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv])
      }
      if(return_lS){
        lS[levellcv, dimlcv] = 2*sum(log(diag(cS)))
        for(paralcv in 1:CGGP$numpara){
          dlS[levellcv, CGGP$numpara*(dimlcv-1)+paralcv] = -sum(diag(dMatdtheta[[(dimlcv-1)*Q+levellcv]][1:nv,nv*(paralcv-1)+1:nv]))
        }
      }
    }
  }
  
  if(is.matrix(y)){
    numout = dim(y)[2]
    sigma2 = rep(0,numout) # Predictive weight for each measured point
    dsigma2 = matrix(0,nrow=CGGP$numpara*CGGP$d,ncol=numout) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        VVV1=unlist(cholS[gg+CGGP$uo[blocklcv,]])
        VVV2=unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]])
        VVV3=CGGP$gridsizest[blocklcv,]
        for(outdimlcv in 1:numout){
          B0 = y[IS,outdimlcv]
          B = (CGGP$w[blocklcv]/dim(y)[1])*B0
          dB = rcpp_gkronDBS(VVV1,VVV2,B,VVV3)
          dsigma2[,outdimlcv] = dsigma2[,outdimlcv] + as.vector(dB%*%B0)
          sigma2[outdimlcv] = sigma2[outdimlcv]  + sum(B0*B)
        }
      }
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }else{
    sigma2 = 0 # Predictive weight for each measured point
    dsigma2 = rep(0,nrow=CGGP$d) # Predictive weight for each measured point
    gg = (1:CGGP$d-1)*Q
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv])>0.5){
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]];
        B0 = y[IS]
        B = (CGGP$w[blocklcv]/length(y))*B0
        dB = rcpp_gkronDBS(unlist(cholS[gg+CGGP$uo[blocklcv,]]),unlist(dMatdtheta[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
        
        dsigma2 = dsigma2 +t(B0)%*%t(dB)
        sigma2 = sigma2 + t(B0)%*%B
      }
    }
    out <- list(sigma2=sigma2,
                dsigma2=dsigma2)
    if (return_lS) {
      out$lS <- lS
      out$dlS <- dlS
    }
  }
  out
}
