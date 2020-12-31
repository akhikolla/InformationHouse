#' Calculate negative log posterior
#'
#' @param theta Correlation parameters
#' @param CGGP CGGP object
#' @param y Measured values of CGGP$design
#' @param ... Forces you to name remaining arguments
#' @param Xs Supplementary input data
#' @param ys Supplementary output data
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' @md
#'
#' @return Likelihood
#' @export
#' @useDynLib CGGP
#'
#' @examples
#' cg <- CGGPcreate(d=3, batchsize=20)
#' Y <- apply(cg$design, 1, function(x){x[1]+x[2]^2})
#' cg <- CGGPfit(cg, Y)
#' CGGP_internal_neglogpost(cg$thetaMAP, CGGP=cg, y=cg$y)
CGGP_internal_neglogpost <- function(theta, CGGP, y, ..., ys=NULL, Xs=NULL,
                                     HandlingSuppData = "Correct") {
  # Return Inf if theta is too large
  if (max(theta) >= 1 || min(theta) <= -1) {
    return(Inf)
  }
  
  # ================================
  # Check that inputs are acceptable
  # ================================
  if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore"))) {
    stop(paste("HandlingSuppData in CGGP_internal_neglogpost must be one of",
               "Correct, Only, Ignore"))
  }
  if(!(is.null(ys) || length(ys)==0) && (is.null(y) || length(y)==0)){
    HandlingSuppData = "Only"
  }else if((is.null(ys) || length(ys)==0) && !(is.null(y) || length(y)==0)){
    HandlingSuppData = "Ignore"
  }else if((is.null(ys) || length(ys)==0) && (is.null(y) || length(y)==0)){
    stop(paste("You have given no y or ys to CGGP_internal_neglogpost"))
  }
  
  
  # =======================
  # Calculate neglogpost
  # =======================
  
  if(HandlingSuppData == "Only" || HandlingSuppData == "Correct"){
    Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv],
                       theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       returnlogs=TRUE)
      Sigma_t = Sigma_t+V
    }
    Sigma_t = exp(Sigma_t)
    
    Sigma_t = (1-CGGP$nugget)*Sigma_t+diag(dim(Sigma_t)[1])*CGGP$nugget
  }
  
  if(HandlingSuppData == "Only"){
    try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
    if (inherits(try.chol, "try-error")) {
      return(Inf)
    }; rm(try.chol)
    
    Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose=TRUE))
    sigma2_hat_supp = colSums(as.matrix(ys*Sti_resid))/dim(Xs)[1]
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
  }
  
  if(HandlingSuppData == "Ignore"){
    sigma2anddsigma2 <- CGGP_internal_calcsigma2(CGGP=CGGP, y=y, theta=theta,
                                                 return_lS=TRUE)
    lS <- sigma2anddsigma2$lS
    sigma2_hat_grid = sigma2anddsigma2$sigma2
    
    lDet_grid = 0
    for (blocklcv in 1:CGGP$uoCOUNT) {
      nv = CGGP$gridsize[blocklcv]/CGGP$gridsizes[blocklcv,]
      uonow = CGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        lDet_grid = lDet_grid + (
          lS[uonow[dimlcv], dimlcv] - lS[uonow[dimlcv] - 1, dimlcv]
        )*nv[dimlcv]
      }
    }
  }
  
  if(HandlingSuppData == "Correct"){
    lik_stuff <- CGGP_internal_calc_cholS_lS_sigma2_pw(CGGP=CGGP, y=y,
                                                       theta=theta)
    cholS = lik_stuff$cholS
    lS <- lik_stuff$lS
    sigma2_hat_grid = lik_stuff$sigma2
    pw = lik_stuff$pw
    
    lDet_grid = 0 # Not needed for glik, only for lik
    
    for (blocklcv in 1:CGGP$uoCOUNT) {
      nv = CGGP$gridsize[blocklcv]/CGGP$gridsizes[blocklcv,]
      uonow = CGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] -
                                   lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
      }
    }
    
    #For these three, I need Cs,pw,allChols
    Cs = (matrix(0,dim(Xs)[1],CGGP$ss))
    GGGG = list(matrix(1,dim(Xs)[1],length(CGGP$xb)),CGGP$d)
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(Xs[,dimlcv],
                       CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]],
                       theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       returnlogs=TRUE)
      GGGG[[dimlcv]] = exp(V)
      Cs = Cs+V[,CGGP$designindex[,dimlcv]]
    }
    Cs = exp(Cs)
    yhats = Cs%*%pw
  }
  
  if(HandlingSuppData == "Correct" ){
    Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv],
                       theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       returnlogs=TRUE)
      Sigma_t = Sigma_t+V
    }
    Sigma_t = exp(Sigma_t)
    
    MSE_s = matrix(NaN,nrow=dim(Xs)[1]*dim(Xs)[1],
                   ncol=(CGGP$d)*(CGGP$maxlevel))
    Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
    for (dimlcv in 1:CGGP$d) {
      gg = (dimlcv-1)*Q
      TT1 = GGGG[[dimlcv]]
      for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
        INDSN = 1:CGGP$sizest[levellcv]
        INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],
                           index.return = TRUE)$ix]
        REEALL = CGGP_internal_postvarmatcalc_fromGMat(
          TT1,
          c(),
          as.matrix(cholS[[gg+levellcv]]),
          c(),
          INDSN,
          CGGP$numpara)
        MSE_s[,(dimlcv-1)*CGGP$maxlevel+levellcv]  =  as.vector(REEALL)
      }
    }
    
    Sigma_t2 = as.vector(Sigma_t)
    rcpp_fastmatclcr(CGGP$uo[1:CGGP$uoCOUNT,], CGGP$w[1:CGGP$uoCOUNT],
                     MSE_s,Sigma_t2,CGGP$maxlevel)
    Sigma_t = matrix(Sigma_t2,nrow=dim(Xs)[1] , byrow = FALSE)
    
    Sigma_t = (1-CGGP$nugget)*Sigma_t+diag(dim(Sigma_t)[1])*CGGP$nugget
    
    try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
    if (inherits(try.chol, "try-error")) {
      return(Inf)
    }; rm(try.chol)
    
    Sti_resid = backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,
                                               transpose = TRUE))
    sigma2_hat_supp = colSums((ys-yhats)*Sti_resid)/dim(Xs)[1]
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
  }
  
  if(HandlingSuppData == "Ignore"){
    sigma2_hat = sigma2_hat_grid
    lDet = lDet_grid
    
    if(!is.matrix(y)){
      nsamples = length(y)
    }else{
      nsamples = dim(y)[1]
    }
  } 
  
  if(HandlingSuppData == "Only"){
    sigma2_hat = sigma2_hat_supp
    lDet = lDet_supp
    
    if(!is.matrix(ys)){
      nsamples = length(ys)
    }else{
      nsamples = dim(ys)[1]
    }
  }
  
  if(HandlingSuppData =="Correct"){
    sigma2_hat = sigma2_hat_grid*dim(CGGP$design)[1] / (
      dim(Xs)[1]+dim(CGGP$design)[1])+sigma2_hat_supp*dim(Xs)[1]/(
        dim(Xs)[1]+dim(CGGP$design)[1])
    lDet = lDet_grid+lDet_supp
    
    if(!is.matrix(y)){
      nsamples = length(y)+length(ys)
    }else{
      nsamples = dim(y)[1]+dim(ys)[1]
    }
  }
  
  neglogpost =  0.1*sum((log(1-theta)-log(theta+1))^2) #start out with prior
  neglogpost =  neglogpost+0.1*sum((log(1-theta)-log(theta+1))^4) #start out with prior
  
  output_is_1D <- if (!is.null(y)) {!is.matrix(y)} else {!is.matrix(ys)}
  if(output_is_1D){
    neglogpost = neglogpost+1/2*(nsamples*log(sigma2_hat[1])+lDet)
  }else{
    neglogpost = neglogpost+1/2*(nsamples*mean(log(c(sigma2_hat)))+lDet)
  }
  
  return(neglogpost)
}
