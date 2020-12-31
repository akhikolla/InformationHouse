

#' Gradient of negative log likelihood posterior
#'
#' @param theta Log of correlation parameters
#' @param CGGP CGGP object
#' @param y CGGP$design measured values
#' @param return_lik If yes, it returns a list with lik and glik
#' @param ... Forces you to name remaining arguments
#' @param Xs Supplementary input data
#' @param ys Supplementary output data
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#'
#' @return Vector for gradient of likelihood w.r.t. x (theta)
#' @export
#'
#' @examples
#' cg <- CGGPcreate(d=3, batchsize=20)
#' Y <- apply(cg$design, 1, function(x){x[1]+x[2]^2})
#' cg <- CGGPfit(cg, Y)
#' CGGP_internal_gneglogpost(cg$thetaMAP, CGGP=cg, y=cg$y)
CGGP_internal_gneglogpost <- function(theta, CGGP, y,..., return_lik=FALSE,
                                      ys=NULL,Xs=NULL,
                                      HandlingSuppData = "Correct") {
  # ================================
  # Check that inputs are acceptable
  # ================================
  if (!(HandlingSuppData %in% c("Correct", "Only", "Ignore"))) {
    stop(paste("HandlingSuppData in CGGP_internal_neglogpost must be one of",
               "Correct, Only, Ignore"))
  }
  
  if(!is.null(Xs) && (is.null(y) || length(y)==0)){
    HandlingSuppData = "Only"
  }else if(is.null(Xs) && !(is.null(y) || length(y)==0)){
    HandlingSuppData = "Ignore"
  }else if(is.null(Xs) && (is.null(y) || length(y)==0)){
    stop(paste("You have given no y or ys to CGGP_internal_neglogpost"))
  }
  
  
  # =======================
  # Calculate gneglogpost
  # =======================
  
  if(HandlingSuppData == "Only" || HandlingSuppData == "Correct"){
    Sigma_t = matrix(0,dim(Xs)[1],dim(Xs)[1])
    dSigma_to = list(Sigma_t,CGGP$d) 
    RTR = list(Sigma_t,CGGP$d) 
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(Xs[,dimlcv], Xs[,dimlcv],
                       theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       return_dCdtheta=TRUE,returnlogs=TRUE)
      dSigma_to[[dimlcv]] = V$dCdtheta
      Sigma_t = Sigma_t+V$C
    }
    Sigma_t = exp(Sigma_t)
    
    Cmat1 = matrix(rep(Sigma_t, CGGP$numpara), nrow = nrow(Sigma_t),
                   byrow = FALSE )
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      dSigma_to[[dimlcv]] =Cmat1*dSigma_to[[dimlcv]]
    }
    
    Sigma_t = (1-CGGP$nugget)*Sigma_t+diag(dim(Sigma_t)[1])*CGGP$nugget
    for (dimlcv in 1:CGGP$d) {
      dSigma_to[[dimlcv]] = (1-CGGP$nugget)*dSigma_to[[dimlcv]]
    }
  }
  
  
  if(HandlingSuppData == "Only"){
    try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
    if (inherits(try.chol, "try-error")) {
      stop("chol error in gneglogpost #1, this can happen when neglogpost is Inf")
      # This came up a lot when running nlminb on the initial point.
      # If the initial neglogpost is Inf, it will call gneglogpost
      # and get the error here. To avoid this we have to make sure the
      # initial points of nlminb are always finite values.
      # This one wasn't the problem, it was usually the other which
      # happens when there is supp data.
      # return(rep(NA, length(theta)))
    }; rm(try.chol)
    
    tempvec1 = backsolve(Sigma_chol,backsolve(Sigma_chol,ys,transpose = TRUE))
    sigma2_hat_supp = colSums(as.matrix(ys*tempvec1))/dim(Xs)[1]
    if(is.matrix(ys)){
      dsigma2_hat_supp = matrix(0,CGGP$d*CGGP$numpara,ncol=dim(ys)[2])
    }else{
      dsigma2_hat_supp = rep(0,CGGP$d*CGGP$numpara)
    }
    
    dlDet_supp=rep(0,CGGP$d*CGGP$numpara)
    for (dimlcv in 1:CGGP$d) {
      for(paralcv in 1:CGGP$numpara){
        dSigma_supp = as.matrix((
          dSigma_to[[dimlcv]])[ ,((paralcv-1)*dim(Sigma_chol)[2]+1
          ):(paralcv*dim(Sigma_chol)[2])])
        tempvec2= dSigma_supp%*%tempvec1
        if(is.matrix(dsigma2_hat_supp )){
          if(dim(dsigma2_hat_supp)[1]>1.5){
            dsigma2_hat_supp[(dimlcv-1)*CGGP$numpara+paralcv,] =
              -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
          }else{
            dsigma2_hat_supp[,(dimlcv-1)*CGGP$numpara+paralcv] =
              -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
          }
        }else{
          dsigma2_hat_supp[(dimlcv-1)*CGGP$numpara+paralcv] =
            -colSums(as.matrix(tempvec1*tempvec2))/dim(Xs)[1]
        }
        dlDet_supp[(dimlcv-1)*CGGP$numpara+paralcv] =
          sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_supp,
                                                  transpose = TRUE))))
      }
    }
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
    
  }
  
  if(HandlingSuppData == "Ignore" ){
    
    sigma2anddsigma2 <- CGGP_internal_calcsigma2anddsigma2(CGGP=CGGP, y=y,
                                                           theta=theta,
                                                           return_lS=TRUE)
    lS <- sigma2anddsigma2$lS
    dlS <-sigma2anddsigma2$dlS
    
    sigma2_hat_grid = sigma2anddsigma2$sigma2
    dsigma2_hat_grid = sigma2anddsigma2$dsigma2
    
    lDet_grid = 0
    dlDet_grid = rep(0, CGGP$numpara*CGGP$d) 
    
    for (blocklcv in 1:CGGP$uoCOUNT) {
      nv = CGGP$gridsize[blocklcv]/CGGP$gridsizes[blocklcv,]
      uonow = CGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        if (return_lik) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] -
                                     lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
        IS = (dimlcv-1)*CGGP$numpara+1:CGGP$numpara
        dlDet_grid[IS] = dlDet_grid[IS] + (dlS[uonow[dimlcv], IS] -
                                             dlS[uonow[dimlcv]-1, IS]
        )*nv[dimlcv]
      }
    }
  }
  if(HandlingSuppData == "Correct"){
    Cs = matrix(0,dim(Xs)[1],CGGP$ss)
    dCs = list(matrix(0,dim(Xs)[1],CGGP$numpara*CGGP$ss),dim(Xs)[2])
    GGGG = list(matrix(1,dim(Xs)[1],length(CGGP$xb)),CGGP$d)
    dGGGG1 = list(matrix(1,dim(Xs)[1],length(CGGP$xb)),CGGP$d)
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      Xbn = CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]]
      V = CGGP$CorrMat(Xs[,dimlcv],  Xbn,
                       theta[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       return_dCdtheta=TRUE,returnlogs=TRUE)
      A = round(as.vector(outer(CGGP$designindex[,dimlcv],
                                length(Xbn)*(0:(CGGP$numpara-1)),
                                FUN=function(x,y) x+y)))
      GGGG[[dimlcv]] = exp(V$C)
      dGGGG1[[dimlcv]] = V$dCdtheta
      dCs[[dimlcv]] = V$dCdtheta[,A]
      Cs= Cs+V$C[,CGGP$designindex[,dimlcv]]
    }
    Cs = exp(Cs)
    
    Cmat1 = matrix(rep(Cs, CGGP$numpara), nrow = nrow(Cs), byrow = FALSE)
    
    dCs = lapply(dCs, FUN = function(x) Cmat1*x)
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      Cmat2 = matrix(rep(GGGG[[dimlcv]], CGGP$numpara),
                     nrow = nrow(GGGG[[dimlcv]]) , byrow = FALSE)
      dGGGG1[[dimlcv]] =Cmat2*(dGGGG1[[dimlcv]])
    }
    
    lik_stuff <- CGGP_internal_calc_cholS_lS_dsigma2_pw_dMatdtheta(
      CGGP=CGGP, y=y, theta=theta)
    cholS = lik_stuff$cholS
    dSV = lik_stuff$dMatdtheta
    lS <- lik_stuff$lS
    dlS <- lik_stuff$dlS
    sigma2_hat_grid = lik_stuff$sigma2
    dsigma2_hat_grid = lik_stuff$dsigma2
    pw = lik_stuff$pw
    
    yhats = Cs%*%pw
    
    lDet_grid = 0
    dlDet_grid = rep(0, CGGP$numpara*CGGP$d) 
    
    for (blocklcv in 1:CGGP$uoCOUNT) {
      nv = CGGP$gridsize[blocklcv]/CGGP$gridsizes[blocklcv,]
      uonow = CGGP$uo[blocklcv,]
      for (dimlcv in which(uonow>1.5)) {
        if (return_lik) {
          lDet_grid = lDet_grid + (lS[uonow[dimlcv], dimlcv] -
                                     lS[uonow[dimlcv] - 1, dimlcv])*nv[dimlcv]
        }
        IS = (dimlcv-1)*CGGP$numpara+1:CGGP$numpara
        dlDet_grid[IS] = dlDet_grid[IS] + (dlS[uonow[dimlcv], IS] -
                                             dlS[uonow[dimlcv]-1, IS]
        )*nv[dimlcv]
      }
    }
  }
  
  
  if(HandlingSuppData == "Correct"){
    MSE_s = matrix(NaN,nrow=dim(Xs)[1]*dim(Xs)[1],
                   ncol=(CGGP$d)*(CGGP$maxlevel))
    dMSE_s = matrix(NaN,nrow=dim(Xs)[1]*dim(Xs)[1],
                    ncol=CGGP$numpara*CGGP$d*CGGP$maxlevel)
    Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
    for (dimlcv in 1:CGGP$d) {
      gg = (dimlcv-1)*Q
      TT1 = GGGG[[dimlcv]]
      TT2  = dGGGG1[[dimlcv]]
      for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
        INDSN = 1:CGGP$sizest[levellcv]
        INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],
                           index.return = TRUE)$ix]
        REEALL = CGGP_internal_postvarmatcalc_fromGMat(TT1,
                                                       TT2,
                                                       cholS[[gg+levellcv]],
                                                       dSV[[gg+levellcv]],
                                                       INDSN,
                                                       CGGP$numpara,
                                                       returndG = TRUE,
                                                       returnderiratio =TRUE)
        MSE_s[,(dimlcv-1)*CGGP$maxlevel+levellcv] = as.vector(REEALL$Sigma_mat)
        for(paralcv in 1:CGGP$numpara){
          dMSE_s[ , CGGP$numpara*(dimlcv-1)*CGGP$maxlevel +
                    (paralcv-1)*CGGP$maxlevel+levellcv] =
            as.vector(REEALL$dSigma_mat[,((dim(Xs)[1]*(
              paralcv-1)+1):(dim(Xs)[1]*paralcv))])
        }
      }
      
    }
    
    dsigma2_hat_part1 = 0*dsigma2_hat_grid
    dsigma2_hat_part2 = 0*dsigma2_hat_grid
    dsigma2_hat_part3 = 0*dsigma2_hat_grid
    dlDet_supp = 0*dlDet_grid
    Sigma_t2 = as.vector(Sigma_t)
    dSigma_to2 = matrix(0,length(Sigma_t2),CGGP$numpara*CGGP$d)
    
    for (dimlcv in 1:CGGP$d) {
      for(paralcv in 1:CGGP$numpara){
        dSigma_to2[,CGGP$numpara*(dimlcv-1)+paralcv] =
          as.vector(t((dSigma_to[[dimlcv]])[,((dim(Xs)[1]*(paralcv-1)+1
          ):(dim(Xs)[1]*paralcv))]))
      }
    }
    
    rcpp_fastmatclcranddclcr(CGGP$uo[1:CGGP$uoCOUNT,], CGGP$w[1:CGGP$uoCOUNT],
                             MSE_s, dMSE_s, Sigma_t2, dSigma_to2,
                             CGGP$maxlevel, CGGP$numpara)
    
    for (dimlcv in 1:CGGP$d) {
      dSigma_to[[dimlcv]] = matrix(dSigma_to2[,CGGP$numpara*(dimlcv-1) +
                                                (1:CGGP$numpara)],
                                   nrow=dim(Xs)[1] , byrow = FALSE)
    }
    Sigma_t = matrix(Sigma_t2,nrow=dim(Xs)[1] , byrow = FALSE)
    
    Sigma_t = (1-CGGP$nugget)*Sigma_t+diag(dim(Sigma_t)[1])*CGGP$nugget
    for (dimlcv in 1:CGGP$d) {
      dSigma_to[[dimlcv]] = (1-CGGP$nugget)*dSigma_to[[dimlcv]]
    }
    
    try.chol <- try({Sigma_chol = chol(Sigma_t)}, silent = TRUE)
    if (inherits(try.chol, "try-error")) {
      # stop(paste("chol error in gneglogpost #2, this can happen when",
      #            " neglogpost is Inf, theta is ", theta, collapse=' '))
      warning(paste(c("chol error in gneglogpost #2, this can happen when",
                 " neglogpost is Inf, theta is ", theta)))
      return(NaN * theta)
      # This came up a lot when running nlminb on the initial point.
      # If the initial neglogpost is Inf, it will call gneglogpost
      # and get the error here. To avoid this we have to make sure the
      # initial points of nlminb are always finite values.
      # return(rep(NA, length(theta)))
    }; rm(try.chol)
    
    tempvec1= backsolve(Sigma_chol,backsolve(Sigma_chol,ys-yhats,
                                             transpose = TRUE))
    
    for (dimlcv in 1:CGGP$d) {
      for(paralcv in 1:CGGP$numpara){
        dCpn = as.matrix((dCs[[dimlcv]])[,((paralcv-1)*dim(Cs)[2]+1
        ):(paralcv*dim(Cs)[2])])
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
            dsigma2_hat_part2[(dimlcv-1)*CGGP$numpara+paralcv,] =
              -2*colSums((tempvec1)*(dCpn%*%pw))
          }else{
            dsigma2_hat_part2[,(dimlcv-1)*CGGP$numpara+paralcv] =
              -2*colSums((tempvec1)*(dCpn%*%pw))
          }
        }else{
          dsigma2_hat_part2[(dimlcv-1)*CGGP$numpara+paralcv] =
            -2*colSums((tempvec1)*(dCpn%*%pw))
        }
        
        
        dSigma_now = as.matrix((
          dSigma_to[[dimlcv]]
        )[,((paralcv-1)*dim(Sigma_chol)[2]+1):(paralcv*dim(Sigma_chol)[2])])
        tempvec2= dSigma_now%*%tempvec1
        if(is.matrix(dsigma2_hat_part2)){
          if(dim(dsigma2_hat_part2)[1]>1.5){
            dsigma2_hat_part1[(dimlcv-1)*CGGP$numpara+paralcv,] =
              -colSums(tempvec1*tempvec2)
          }else{
            dsigma2_hat_part1[,(dimlcv-1)*CGGP$numpara+paralcv] =
              -colSums(tempvec1*tempvec2)
          }
        }else{
          dsigma2_hat_part1[(dimlcv-1)*CGGP$numpara+paralcv] =
            -colSums(tempvec1*tempvec2)
        }
        
        
        dlDet_supp[(dimlcv-1)*CGGP$numpara+paralcv] =
          sum(diag(backsolve(Sigma_chol,backsolve(Sigma_chol,dSigma_now,
                                                  transpose = TRUE))))
      }
    }
    
    if(is.vector(y)){
      temp4 = as.vector(t(Cs)%*%tempvec1)
    }else{
      temp4 = t(Cs)%*%tempvec1
    }
    dsigma2_hat_part3 =  -2*(CGGP_internal_calc_dvalo(CGGP,y,temp4,cholS,dSV
    )$dvalo)
    lDet_supp = 2*sum(log(diag(Sigma_chol)))
    sigma2_hat_supp = colSums((ys-yhats)*tempvec1)/dim(Xs)[1]
    dsigma2_hat_supp = (dsigma2_hat_part1+dsigma2_hat_part2 +
                          dsigma2_hat_part3)/dim(Xs)[1]
  }
  
  
  
  if(HandlingSuppData == "Ignore"){
    sigma2_hat = sigma2_hat_grid
    dsigma2_hat = dsigma2_hat_grid
    dlDet = dlDet_grid
    lDet = lDet_grid
    
    if(!is.matrix(y)){
      nsamples = length(y)
    }else{
      nsamples = dim(y)[1]
      ndim = dim(y)[2]
    }
  } 
  
  if(HandlingSuppData == "Only"){
    sigma2_hat = sigma2_hat_supp
    dsigma2_hat = dsigma2_hat_supp
    dlDet = dlDet_supp
    lDet = lDet_supp
    
    if(!is.matrix(ys)){
      nsamples = length(ys)
    }else{
      nsamples = dim(ys)[1]
      ndim = dim(ys)[2]
    }
  }
  
  if(HandlingSuppData =="Correct"){
    sigma2_hat = sigma2_hat_grid *
      dim(CGGP$design)[1]/(dim(Xs)[1] +dim(CGGP$design)[1]) +
      sigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(CGGP$design)[1])
    dsigma2_hat = dsigma2_hat_grid *
      dim(CGGP$design)[1]/(dim(Xs)[1]+dim(CGGP$design)[1]) +
      dsigma2_hat_supp*dim(Xs)[1]/(dim(Xs)[1]+dim(CGGP$design)[1])
    dlDet = dlDet_grid+dlDet_supp
    lDet = lDet_grid+lDet_supp
    
    if(!is.matrix(y)){
      nsamples = length(y)+length(ys)
    }else{
      nsamples = dim(y)[1]+dim(ys)[1]
      ndim = dim(y)[2]
    }
  }
  neglogpost =  0.1*sum((log(1-theta)-log(theta+1))^2) #start out with prior
  gneglogpost = -0.2*(log(1-theta)-log(theta+1))*((1/(1-theta))+1/(1+theta))
  
  neglogpost =  neglogpost+0.1*sum((log(1-theta)-log(theta+1))^4) #start out with prior
  gneglogpost = gneglogpost-0.1*4*((log(1-theta)-log(theta+1))^3)*((1/(1-theta))+1/(1+theta))
  
  output_is_1D <- if (!is.null(y)) {!is.matrix(y)} else {!is.matrix(ys)}
  if(output_is_1D){
    neglogpost =neglogpost +1/2*(nsamples*log(sigma2_hat[1])+lDet)
    gneglogpost = gneglogpost+1/2*(dlDet+nsamples*dsigma2_hat / sigma2_hat[1])
  }else{
    neglogpost = neglogpost+1/2*(nsamples*mean(log(c(sigma2_hat)))+lDet)
    gneglogpost = gneglogpost+1/2*dlDet
    for(i in 1:ndim){
      gneglogpost = gneglogpost+1/2*1/ndim*nsamples*dsigma2_hat[,i]/sigma2_hat[i]
    }
  }
  
  
  if(return_lik){
    return(list(neglogpost=neglogpost,gneglogpost=gneglogpost))
  } else {
    return(gneglogpost)
  }
}
