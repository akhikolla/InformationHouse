#' Predict with CGGP object
#' 
#' Predict using SG with y values at xp?
#' Shouldn't y values already be stored in SG?
#'
#' @param xp x value to predict at
#' @param CGGP SG object
#' @param theta Leave as NULL unless you want to use a value other than thetaMAP.
#' Much slower.
#' @param outdims If multiple outputs fit without PCA and with separate
#' parameters, you can predict just for certain dimensions to speed it up.
#' Will leave other columns in the output, but they will be wrong.
#'
#' @return Predicted mean values
#' @export
#' @family CGGP core functions
#'
#' @examples
#' SG <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
#' SG <- CGGPfit(SG, Y=y)
#' CGGPpred(SG, matrix(c(.1,.1,.1),1,3))
#' cbind(CGGPpred(SG, SG$design)$mean, y) # Should be near equal
CGGPpred <- function(CGGP, xp, theta=NULL, outdims=NULL) {
  if (!inherits(CGGP, "CGGP")) {
    stop("First argument to CGGP must be an CGGP object")
  }
  # Require that you run CGGPfit first
  if (is.null(CGGP$supplemented)) {
    stop("You must run CGGPfit on CGGP object before using CGGPpredict")
  }
  if (CGGP$supplemented && (is.null(CGGP[["Y"]]) || length(CGGP$Y)==0)) {
    return(CGGP_internal_predwithonlysupp(CGGP=CGGP, xp=xp, theta=theta, outdims=outdims))
  }
  # We could check for design_unevaluated, maybe give warning?
  
  # If xp has many rows, split it up over calls to CGGPpred and regroup
  predgroupsize <- 100
  if (nrow(xp) >= predgroupsize*2) {
    ngroups <- floor(nrow(xp) / predgroupsize)
    list_preds <- lapply(1:ngroups,
                         function(i) {
                           inds.i <- if (i < ngroups) {1:predgroupsize + (i-1)*predgroupsize} else {((i-1)*predgroupsize+1):(nrow(xp))}
                           CGGPpred(CGGP=CGGP,
                                    xp=xp[inds.i, , drop=FALSE],
                                    theta=theta,
                                    outdims=outdims)
                         })
    outlist <- list(mean=do.call(rbind, lapply(list_preds, function(ll) ll$mean)),
                    var= do.call(rbind, lapply(list_preds, function(ll) ll$var ))
    )
    if (nrow(outlist$mean)!=nrow(xp) || nrow(outlist$var)!=nrow(xp)) {
      stop("Error with CGGPpred predicting many points")
    }
    return(outlist)
  }
  
  # If theta is given (for full Bayesian prediction), need to recalculate pw
  if (!is.null(theta) && length(theta)!=length(CGGP$thetaMAP)) {stop("Theta is wrong length")}
  if (!is.null(theta) && all(theta==CGGP$thetaMAP)) {
    # If you give in theta=thetaMAP, set it to NULL to avoid recalculating.
    theta <- NULL
  }
  if (is.null(theta)) {
    thetaMAP <- CGGP$thetaMAP
    recalculate_pw <- FALSE
  } else {
    thetaMAP <- theta
    rm(theta)
    # pw <- CGGP_internal_calcpw(CGGP, CGGP$y, theta=thetaMAP)
    recalculate_pw <- TRUE
  }
  
  
  separateoutputparameterdimensions <- is.matrix(CGGP$thetaMAP)
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions) {
    ncol(CGGP$y)
  } else {
    1
  }
  if (nopd > 1) {
    # meanall <- matrix(NaN, nrow(xp), ncol=nopd)
    meanall2 <- matrix(0, nrow(xp), ncol=ncol(CGGP$Y))
    # varall <- matrix(NaN, nrow(xp), ncol=ncol(CGGP$Y))
    tempvarall <- matrix(0, nrow(xp), ncol=ncol(CGGP$Y))
  }
  
  if (!is.null(outdims) && nopd==1) {
    stop("outdims can only be given when multiple outputs and separate correlation parameters")
  }
  opd_values <- if (is.null(outdims)) {1:nopd} else {outdims}
  for (opdlcv in opd_values) {# 1:nopd) {
    thetaMAP.thisloop <- if (nopd==1) thetaMAP else thetaMAP[, opdlcv]
    if (!recalculate_pw) { # use already calculated
      pw.thisloop <- if (nopd==1) CGGP$pw else CGGP$pw[,opdlcv]
      sigma2MAP.thisloop <- CGGP$sigma2MAP
      cholS.thisloop <- if (nopd==1) CGGP$cholSs else CGGP$cholSs[[opdlcv]]
    } else { # recalculate pw and sigma2MAP
      y.thisloop <- if (nopd==1) CGGP$y else CGGP$y[,opdlcv]
      
      lik_stuff <- CGGP_internal_calc_cholS_lS_sigma2_pw(CGGP=CGGP,
                                                         y=y.thisloop,
                                                         theta=thetaMAP.thisloop
      )
      cholS.thisloop = lik_stuff$cholS
      sigma2MAP.thisloop <-  lik_stuff$sigma2
      pw.thisloop = lik_stuff$pw
      # It can be vector when there is multiple output
      sigma2MAP.thisloop <- as.vector(sigma2MAP.thisloop)
      rm(y.thisloop)
    }
    mu.thisloop <- if (nopd==1) CGGP$mu else CGGP$mu[opdlcv] # Not used for PCA, added back at end
    
    # Cp is sigma(x_0) in paper, correlation vector between design points and xp
    Cp = matrix(0,dim(xp)[1],CGGP$ss)
    GGGG = list(matrix(1,dim(xp)[1],length(CGGP$xb)),CGGP$d)
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(xp[,dimlcv], CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]],
                       thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                       returnlogs=TRUE)
      GGGG[[dimlcv]] = exp(V)
      Cp = Cp+V[,CGGP$designindex[,dimlcv]]
    }
    Cp = exp(Cp)
    
    
    ME_t = matrix(1,dim(xp)[1],1)
    MSE_v = list(matrix(0,dim(xp)[1],2),(CGGP$d+1)*(CGGP$maxlevel+1)) 
    Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
    for (dimlcv in 1:CGGP$d) {
      for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
        Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
        gg = (dimlcv-1)*Q
        INDSN = 1:CGGP$sizest[levellcv]
        INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],
                           index.return = TRUE)$ix]
        MSE_v[[(dimlcv)*CGGP$maxlevel+levellcv]] =
          CGGP_internal_postvarmatcalc_fromGMat(GGGG[[dimlcv]],
                                                c(),
                                                as.matrix(
                                                  cholS.thisloop[[gg+levellcv]]
                                                ),
                                                c(),
                                                INDSN,
                                                CGGP$numpara,
                                                returndiag=TRUE)
      }
    }
    
    for (blocklcv in 1:CGGP$uoCOUNT) {
      if(abs(CGGP$w[blocklcv]) > 0.5){
        ME_s = matrix(1,nrow=dim(xp)[1],1)
        for (dimlcv in 1:CGGP$d) {
          levelnow = CGGP$uo[blocklcv,dimlcv]
          ME_s = ME_s*MSE_v[[(dimlcv)*CGGP$maxlevel+levelnow]]
        }
        ME_t = ME_t-CGGP$w[blocklcv]*ME_s
      }
    }
    
    
    if (!CGGP$supplemented) {  
      
      # Return list with mean and var predictions
      if(is.vector(pw.thisloop)){
        if (nopd == 1) {
          mean = (mu.thisloop+Cp%*%pw.thisloop)
          var=sigma2MAP.thisloop*ME_t
        }
        
        # With sepparout and PCA (or not), do this
        if (nopd > 1) {
          # meanall2 <- meanall2 + outer(c(Cp%*%pw.thisloop), CGGP$M[opdlcv, ])
          meanall2[,opdlcv] <- c(Cp%*%pw.thisloop)
          
          # This should be correct variance. Needs to be tested better.
          # Pick out the current dimension, set other values to zero
          tempM <- diag(nopd) #CGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          tempvar <- (as.vector(ME_t)%*%t(diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
        }
      }else{ # y was a matrix
        if(length(sigma2MAP.thisloop)==1){
          stop("When is it a matrix but sigma2MAP a scalar???")
        }else{
          mean = (matrix(rep(mu.thisloop,each=dim(xp)[1]), ncol=ncol(CGGP$Y), byrow=FALSE) +
                    (Cp%*%pw.thisloop))
          var=as.vector(ME_t)%*%t(sigma2MAP.thisloop)
        }
      }
    } else { # CGGP$supplemented is TRUE
      if (!recalculate_pw) {
        pw_uppad.thisloop <- if (nopd==1) CGGP$pw_uppad else CGGP$pw_uppad[,opdlcv]
        supppw.thisloop <- if (nopd==1) CGGP$supppw else CGGP$supppw[,opdlcv]
        Sti.thisloop <- if (nopd==1) CGGP$Sti else CGGP$Sti[,,opdlcv]
      } else {
        stop("Give theta in not implemented in CGGPpred. Need to fix sigma2MAP here too!")
      }
      
      Cps = matrix(0,dim(xp)[1],dim(CGGP$Xs)[1])
      GGGG2 = list(matrix(0,nrow=dim(xp)[1],ncol=dim(CGGP$Xs)[1]),CGGP$d)
      for (dimlcv in 1:CGGP$d) { # Loop over dimensions
        V = CGGP$CorrMat(xp[,dimlcv], 
                         CGGP$Xs[,dimlcv], 
                         thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                         returnlogs=TRUE)
        Cps = Cps+V
        
        V = CGGP$CorrMat(CGGP$Xs[,dimlcv], 
                         CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]], 
                         thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                         returnlogs=TRUE)
        GGGG2[[dimlcv]] = exp(V)
      }
      Cps = exp(Cps)
      
      yhatp = Cp%*%pw_uppad.thisloop + Cps%*%supppw.thisloop
      
      
      MSE_ps = matrix(NaN,nrow=dim(CGGP$Xs)[1]*dim(xp)[1],ncol=(CGGP$d)*(CGGP$maxlevel))
      Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
      for (dimlcv in 1:CGGP$d) {
        gg = (dimlcv-1)*Q
        for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
          INDSN = 1:CGGP$sizest[levellcv]
          INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],index.return = TRUE)$ix]
          REEALL= CGGP_internal_postvarmatcalc_fromGMat_asym(GGGG[[dimlcv]],
                                                             GGGG2[[dimlcv]],
                                                             as.matrix(cholS.thisloop[[gg+levellcv]]),
                                                             INDSN)
          MSE_ps[,(dimlcv-1)*CGGP$maxlevel+levellcv]  =  as.vector(REEALL)
        }
      }
      
      
      Cps2 = as.vector(Cps)
      rcpp_fastmatclcr(CGGP$uo[1:CGGP$uoCOUNT,], CGGP$w[1:CGGP$uoCOUNT], MSE_ps, Cps2,CGGP$maxlevel)
      Cps = matrix(Cps2,ncol=dim(CGGP$Xs)[1] , byrow = FALSE)
      
      
      ME_adj = rowSums((Cps%*%Sti.thisloop)*Cps)
      
      
      ME_t = ME_t-ME_adj
      # Return list with mean and var predictions
      if(is.vector(pw.thisloop)){
        if (nopd == 1) {
          mean = (CGGP$mu+ yhatp)
          if (length(CGGP$sigma2MAP)>1) {warning("If this happens, you should fix var here")}
          # Should this be sigma2MAP.thisloop?
          var=CGGP$sigma2MAP[1]*ME_t
        }
        
        # With sepparout and PCA (or not), do this
        if (nopd > 1) {
          # meanall2 <- meanall2 + outer(c(yhatp), CGGP$M[opdlcv,])
          meanall2[,opdlcv] <- yhatp
          leftvar <- if (is.null(CGGP$leftover_variance)) {0} else {CGGP$leftover_variance}
          tempM <- diag(nopd) #CGGP$M
          tempM[-opdlcv,] <- 0
          tempsigma2.thisloop <- sigma2MAP.thisloop
          tempsigma2.thisloop[-opdlcv] <- 0
          tempvar <- (as.vector(ME_t)%*%t(
            leftvar + diag(t(tempM)%*%diag(tempsigma2.thisloop)%*%(tempM))))
        }
        
      }else{
        if(length(CGGP$sigma2MAP)==1){
          stop("This should never happen #952570")
        }else{
          leftvar <- if (is.null(CGGP$leftover_variance)) {0} else {CGGP$leftover_variance}
          mean = matrix(rep(CGGP$mu,each=dim(xp)[1]), ncol=ncol(CGGP$Y), byrow=FALSE)+ yhatp
          var=as.vector(ME_t)%*%t(leftvar+CGGP$sigma2MAP)
        }
      }
      
      
    }
    rm(Cp,ME_t, MSE_v, V) # Just to make sure nothing is carrying through
    if (nopd > 1) {tempvarall <- tempvarall + tempvar}
  }
  # If PCA values were calculated separately, need to do transformation on
  #  both before mu is added, then add mu back
  # if (nopd > 1) {meanall <- sweep(sweep(meanall,2,CGGP$mu) %*% CGGP$M,2,CGGP$mu, `+`)}
  if (nopd > 1) {meanall2 <- sweep(meanall2, 2, CGGP$mu, `+`)}
  
  if (nopd > 1) {
    GP <- list(mean=meanall2, var=tempvarall)
  } else {
    GP <- list(mean=mean, var=var)
  }
  # Check for negative variances, set them all to be a tiny number
  GP$var <- pmax(GP$var, .Machine$double.eps)
  return(GP)
}

#' Calculate MSE prediction along a single dimension
#'
#' @param xp Points at which to calculate MSE
#' @param xl Levels along dimension, vector???
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE predictions
#' @export
#'
#' @examples
#' CGGP_internal_MSEpredcalc(c(.4,.52), c(0,.25,.5,.75,1), theta=c(.1,.2),
#'              CorrMat=CGGP_internal_CorrMatCauchySQ)
CGGP_internal_MSEpredcalc <- function(xp,xl,theta,CorrMat) {
  S = CorrMat(xl, xl, theta)
  n = length(xl)
  cholS = chol(S)
  
  Cp = CorrMat(xp, xl, theta)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  
  MSE_val = 1 - rowSums(t(CiCp)*((Cp)))
  return(MSE_val)
}



#' Calculate posterior variance, faster version
#'
#' @param GMat1 Matrix 1
#' @param GMat2 Matrix 2
#' @param cholS Cholesky factorization of S
#' @param INDSN Indices, maybe
#'
#' @return Variance posterior
## @export
#' @noRd
CGGP_internal_postvarmatcalc_fromGMat_asym <- function(GMat1,GMat2,cholS,INDSN) {
  CoinvC1o = backsolve(cholS,
                       backsolve(cholS,t(GMat1[,INDSN]), transpose = TRUE))
  Sigma_mat = (t(CoinvC1o)%*%(t(GMat2[,INDSN])))
  
  return(Sigma_mat)
  
}
