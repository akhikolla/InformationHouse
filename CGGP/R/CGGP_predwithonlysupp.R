#' Predict with only supplemental data
#'
#' @param CGGP CGGP object with supplemental data
#' @param xp Matrix of points to predict at
#' @param theta Correlation parameters
#' @param outdims Output dimensions to predict for
#'
#' @return Predictions
# @export
#' @noRd
#'
#' @examples
#' d <- 3
#' n <- 30
#' Xs <- matrix(runif(d*n), n, d)
#' Ys <- apply(Xs, 1, function(x){x[1]/(x[3]+.2)+exp(x[3])*cos(x[2]^2)})
#' cg <- CGGPcreate(d, Xs=Xs, Ys=Ys, batchsize=0)
#' xp <- matrix(runif(d*10), ncol=d)
#' predict(cg, xp)
CGGP_internal_predwithonlysupp <- function(CGGP, xp, theta=NULL, outdims=NULL) {
  if (!inherits(CGGP, "CGGP")) {
    stop("First argument to CGGP must be an CGGP object")
  }
  # Require that you run CGGPfit first
  if (is.null(CGGP$supplemented)) {
    stop("You must run CGGPfit on CGGP object before using CGGPpredict")
  }
  if (!CGGP$supplemented) {
    stop("Must be supplemented to use CGGPpred_supponly")
  }
  # We could check for design_unevaluated, maybe give warning?
  
  # If theta is given, need to recalculate pw
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
    ncol(CGGP$ys)
  } else {
    1
  }
  if (nopd > 1) {
    meanall2 <- matrix(0, nrow(xp), ncol=ncol(CGGP$Ys))
    tempvarall <- matrix(0, nrow(xp), ncol=ncol(CGGP$Ys))
  }
  
  if (!is.null(outdims) && nopd==1) {
    stop("outdims can only be given when multiple outputs and separate correlation parameters")
  }
  opd_values <- if (is.null(outdims)) {1:nopd} else {outdims}
  for (opdlcv in opd_values) {# 1:nopd) {
    thetaMAP.thisloop <- if (nopd==1) thetaMAP else thetaMAP[, opdlcv]
    if (!recalculate_pw) { # use already calculated
      sigma2MAP.thisloop <- CGGP$sigma2MAP
      supppw.thisloop <- if (nopd==1) CGGP$supppw else CGGP$supppw[,opdlcv]
      Sti.thisloop <- if (nopd==1) CGGP$Sti else CGGP$Sti[,,opdlcv]
    } else { # recalculate pw and sigma2MAP
      ys.thisloop <- if (nopd==1) CGGP$ys else CGGP$ys[,opdlcv]
      
      Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
      for (dimlcv in 1:CGGP$d) { # Loop over dimensions
        V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
        Sigma_t = Sigma_t*V
      }
      
      Sti_chol <- chol(Sigma_t + diag(CGGP$nugget, nrow(Sigma_t), ncol(Sigma_t)))
      Sti.thisloop <- chol2inv(Sti_chol)
      
      
      supppw.thisloop <- backsolve(Sti_chol, backsolve(Sti_chol, ys.thisloop, transpose = T))
      if (is.matrix(supppw.thisloop) && ncol(supppw.thisloop)==1) {
        supppw.thisloop <- as.vector(supppw.thisloop)
      }
      
      
      sigma2MAP.thisloop <- (t(ys.thisloop) %*% supppw.thisloop) / nrow(CGGP$Xs)
      if (is.matrix(sigma2MAP.thisloop)) {sigma2MAP.thisloop <- diag(sigma2MAP.thisloop)}
      
      rm(Sigma_t, V, Sti_chol, ys.thisloop)
    }
    mu.thisloop <- if (nopd==1) CGGP$mu else CGGP$mu[opdlcv] # Not used for PCA, added back at end
    
    Cps = matrix(1,dim(xp)[1],dim(CGGP$Xs)[1])
    for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      V = CGGP$CorrMat(xp[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
      Cps = Cps*V
    }
    
    yhatp = Cps%*%supppw.thisloop
    
    # Return list with mean and var predictions
    if(is.vector(supppw.thisloop)){
      if (nopd == 1) {
        mean = (CGGP$mu + yhatp)
        var=CGGP$sigma2MAP[1]* (1-diag(Cps %*% CGGP$Sti %*% t(Cps)))  # ME_t
        # could return cov mat here
      }
      
      # With sepparout and PCA (or not), do this
      if (nopd > 1) {
        meanall2[,opdlcv] <- yhatp
        tempvar <- matrix(0, nrow(xp), nopd)
        tempvar[,opdlcv] <- CGGP$sigma2MAP[opdlcv]* (1-diag(Cps %*% CGGP$Sti[,,opdlcv] %*% t(Cps)))
      }
      
    }else{ # supppw is matrix, so predicting multiple columns at once
      if(length(CGGP$sigma2MAP)==1){
        stop("This should never happen #952570")
      }else{
        mean <- matrix(rep(CGGP$mu,each=dim(xp)[1]), ncol=ncol(CGGP$Ys), byrow=FALSE) + yhatp
        var <- sweep(matrix((1-diag(Cps %*% CGGP$Sti %*% t(Cps))), nrow(xp),
                            ncol(CGGP$Ys), byrow = F),
                     2, CGGP$sigma2MAP, `*`)
      }
    }
    
    
    
    # rm(Cp,ME_t, MSE_v, V) # Just to make sure nothing is carrying through
    if (nopd > 1) {tempvarall <- tempvarall + tempvar}
  }
  # If PCA values were calculated separately, need to do transformation on both before mu is added, then add mu back
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
