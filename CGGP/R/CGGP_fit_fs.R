

#' Update CGGP model given data
#' 
#' This function will update the GP parameters for a CGGP design.
#'
#' @param CGGP Sparse grid objects
#' @param Y Output values calculated at CGGP$design
#' @param Xs Supplemental X matrix
#' @param Ys Supplemental Y values
#' @param theta0 Initial theta
#' @param Ynew Values of `CGGP$design_unevaluated`
#' @param separateoutputparameterdimensions If multiple output dimensions,
#' should separate parameters be fit to each dimension?
#' @param set_thetaMAP_to Value for thetaMAP to be set to
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' @param corr Will update correlation function, if left missing it will be
#' same as last time.
#' 
#' @importFrom stats optim rnorm runif nlminb
#'
#' @return Updated CGGP object fit to data given
#' @export
#' @family CGGP core functions
#'
#' @examples
#' cg <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(cg$design, 1, function(x){x[1]+x[2]^2})
#' cg <- CGGPfit(CGGP=cg, Y=y)
CGGPfit <- function(CGGP, Y, Xs=NULL,Ys=NULL,
                    theta0 = pmax(pmin(CGGP$thetaMAP,0.8),-0.8), #gotta pull away from edges to get not corner solution
                    HandlingSuppData=CGGP$HandlingSuppData,
                    separateoutputparameterdimensions=is.matrix(CGGP$thetaMAP),
                    set_thetaMAP_to,
                    corr,
                    Ynew) {
  # ========================================.
  # ==== Check inputs, get Y from Ynew  ====
  # ========================================.
  # If different correlation function is given, update it
  if (!missing(corr)) {
    message("Changing correlation function")
    CGGP <- CGGP_internal_set_corr(CGGP, corr)
  }
  
  # If Y or Ynew is matrix with 1 column, convert it to vector to avoid issues
  if (!missing(Y) && is.matrix(Y) && ncol(Y)==1) {
    Y <- c(Y)
  }
  if (!missing(Ynew) && is.matrix(Ynew) && ncol(Ynew)==1) {
    Ynew <- c(Ynew)
  }
  
  # If Ynew is given, it is only the points that were added last iteration.
  #   Append it to previous Y
  if (!missing(Ynew)) {
    if (!missing(Y)) {stop("Don't give both Y and Ynew, only one")}
    if (is.null(CGGP$Y) || length(CGGP$Y)==0) {
      if (is.matrix(Ynew) && nrow(Ynew) != nrow(CGGP$design_unevaluated)) {
        stop("nrow(Ynew) doesn't match")
      }
      if (!is.matrix(Ynew) && length(Ynew) != nrow(CGGP$design_unevaluated)) {
        stop("length(Ynew) doesn't match")
      }
      Y <- Ynew
    } else if (is.matrix(CGGP$Y)) {
      if (!is.matrix(Ynew)) {stop("Ynew should be a matrix")}
      if (nrow(Ynew) != nrow(CGGP$design_unevaluated)) {
        stop("Ynew is wrong size")
      }
      Y <- rbind(CGGP$Y, Ynew)
    } else { # is numeric vector
      if (length(Ynew) != nrow(CGGP$design_unevaluated)) {
        stop("Ynew is wrong size")
      }
      Y <- c(CGGP$Y, Ynew)
    }
  }
  
  if ((is.matrix(Y) && nrow(Y) == nrow(CGGP$design)) ||
      (length(Y) == nrow(CGGP$design))) {
    CGGP$design_unevaluated <- NULL
  } else {
    stop("CGGP$design and Y have different length")
  }
  
  # ====================================================================.
  # Do the pre-processing                                           ====
  # For cleanness: Y is always the user input, y is after transformation
  # ====================================================================.
  CGGP$Y = Y
  if (any(is.na(Y))) {
    message(paste0(sum(is.na(Y)), "/",length(Y)," Y values are NA, will impute them"))
  }
  if(is.null(Xs)){ # No supplemental data
    CGGP$supplemented = FALSE
    
    if(!is.matrix(Y)){
      CGGP$mu = mean(Y[!is.na(Y)])
      y = Y-CGGP$mu
    }else{ # Y is matrix, PCA no longer an option
      CGGP$mu = colMeans(Y)
      for(oplcv in 1:dim(Y)[2]){
        CGGP$mu[oplcv] = mean(Y[!is.na(Y[,oplcv]),oplcv])
      }
      y <- sweep(Y, 2, CGGP$mu)
      # Need to set CGGP$M somewhere so that it doesn't use transformation
    }
    CGGP$y = y
    
  } else{ # Has supp data, used for prediction but not for fitting params
    # stop("Not working for supp")
    CGGP$supplemented = TRUE
    CGGP$Xs = Xs
    CGGP$Ys = Ys
    
    if(!is.matrix(Y)){
      CGGP$mu = mean(Ys[!is.na(Ys)])
      y = Y-CGGP$mu
      ys = Ys-CGGP$mu
    } else{ # PCA no longer an option
      CGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)),
      for(oplcv in 1:dim(Ys)[2]){
        CGGP$mu[oplcv] = mean(Ys[!is.na(Ys[,oplcv]),oplcv])
      }
      #  or make sure Ys is big enough for this
      y <- sweep(Y, 2, CGGP$mu)
      ys <- sweep(Ys, 2, CGGP$mu)
    }
    CGGP$y = y
    CGGP$ys = ys
  }
  
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions && is.matrix(y)) {
    ncol(y)
  } else {
    1
  }
  
  # Fix theta0
  if (nopd > 1) {
    if (is.vector(theta0)) {
      theta0 <- matrix(theta0, nrow=length(theta0), ncol=nopd, byrow=F)
    }
  }
  
  # Can get an error for theta0 if number of PCA dimensions has changed
  if (is.matrix(theta0) && (ncol(theta0) != nopd)) {
    stop(paste("theta0 should have", nopd, "columns"))
  }
  
  
  # =======================================================.
  # Fit parameters for each output parameter dimension ====
  # =======================================================.
  for (opdlcv in 1:nopd) { # output parameter dimension
    y.thisloop <- if (nopd==1) {y} else {y[,opdlcv]} # All of y or single column
    if (!is.null(Ys)) {ys.thisloop <- if (nopd==1) {ys} else {ys[,opdlcv]}}
    else {ys.thisloop <- NULL}
    theta0.thisloop <- if (nopd==1) {theta0} else {theta0[,opdlcv]}
    
    if(any(is.na(y.thisloop))){
      if (!missing(set_thetaMAP_to) && !is.null(set_thetaMAP_to)) {
        opt.out <- list(par = if (nopd>1) {set_thetaMAP_to[,opdlcv]} else {set_thetaMAP_to})
        y.thisloop=CGGP_internal_imputesomegrid(CGGP,y.thisloop,opt.out)
        repeattimes = 1
      }else{
        y.orig = y.thisloop
        y.thisloop=CGGP_internal_imputesomegrid(CGGP,y.thisloop,theta0.thisloop)
        repeattimes = 10
      }
    }else{
      repeattimes = 1
    }
    
    for(imputelcv in 1:repeattimes){
      if(imputelcv > 1.5){
        if(imputelcv < 2.5){
          y.thisloop=CGGP_internal_imputesomegrid(CGGP,y.orig,thetaMAP)
        }else{
          ystart = y.thisloop
          y.thisloop=CGGP_internal_imputesomegrid(CGGP,y.orig,thetaMAP,ystart=ystart)
          if(max(abs(y.thisloop-ystart))<10^(-10)*max(abs(ystart))){
            break
          }
        }
      }
      # Find MAP theta
      if (!missing(set_thetaMAP_to) && !is.null(set_thetaMAP_to)) {
        opt.out <- list(par = if (nopd>1) {set_thetaMAP_to[,opdlcv]} else {set_thetaMAP_to})
      } else if (is.null(CGGP$Xs)){ # No supp data, just optimize
        opt.out = nlminb(
          theta0.thisloop,
          objective = CGGP_internal_neglogpost,
          gradient = CGGP_internal_gneglogpost,
          y = y.thisloop,
          CGGP = CGGP,
          ys = ys.thisloop,
          Xs = Xs,
          HandlingSuppData=HandlingSuppData,
          control = list(rel.tol = 1e-4,iter.max = 500)
        )
      } else { # W/ supp data, optimize on grid first, then with both
        # Only grid data b/c it's fast
        opt.out = nlminb(
          theta0.thisloop,
          objective = CGGP_internal_neglogpost,
          gradient = CGGP_internal_gneglogpost,
          y = y.thisloop,
          CGGP = CGGP,
          HandlingSuppData="Ignore", # Never supp data here, so set to Ignore
          #  regardless of user setting
          lower=rep(-.9, CGGP$d),
          upper=rep( .9, CGGP$d),
          control = list(rel.tol = 1e-2,iter.max = 500)
        )
        
        neglogpost_par <- CGGP_internal_neglogpost(theta=opt.out$par,
                                                   CGGP=CGGP,
                                                   y=y.thisloop,
                                                   ys=ys.thisloop,
                                                   Xs=Xs,
                                                   HandlingSuppData=HandlingSuppData
        )
        if (is.infinite(neglogpost_par)) {
          theta0_2 <- rep(0, CGGP$d)
        } else {
          theta0_2 <- opt.out$par
        }
        
        # Then use best point as initial point with supp data
        opt.out = nlminb(
          theta0_2,
          objective = CGGP_internal_neglogpost,
          gradient = CGGP_internal_gneglogpost,
          y = y.thisloop,
          ys = ys.thisloop,
          Xs = Xs,
          CGGP = CGGP,
          HandlingSuppData = HandlingSuppData,
          control = list(rel.tol = 1e-4,iter.max = 500)
        )
        
        
      } # End supp data opt
      
      
      thetaMAP <- opt.out$par
      sigma2MAP <- CGGP_internal_calcsigma2anddsigma2(CGGP=CGGP, y=y.thisloop,
                                                      theta=thetaMAP,
                                                      return_lS=FALSE)$sigma2
      
      if(imputelcv > 1.5){
        if(all(abs(sigma2MAP0-sigma2MAP)<0.025*sigma2MAP)){
          break
        }
        sigma2MAP0 = sigma2MAP
      }else{
        # On first time through
        sigma2MAP0 = sigma2MAP
      }
    }
    
    
    # ===================================.
    #  Update parameters and samples ====
    # ===================================.
    
    # Set new theta
    # If one value, it gives it as matrix. Convert it to scalar
    if (length(sigma2MAP) == 1) {sigma2MAP <- sigma2MAP[1,1]}
    
    lik_stuff <- CGGP_internal_calc_cholS_lS_sigma2_pw(CGGP=CGGP, y.thisloop, theta=thetaMAP)
    cholSs = lik_stuff$cholS
    pw <- lik_stuff$pw
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = CGGP_internal_gneglogpost(thetav,CGGP,y.thisloop,
                                      Xs=Xs, ys=ys.thisloop,
                                      HandlingSuppData=HandlingSuppData) *
      (2*(exp(PSTn))/(exp(PSTn)+1)^2)
    
    dimensions_that_need_fixing <- c()
    for(ipara in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[ipara] =10^(-3)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      
      PSTn2=  log((1+thetaMAP)/(1-thetaMAP)) - rsad
      thetav2=(exp(PSTn2)-1)/(exp(PSTn2)+1)
      
      # There can be issues if gneglogpost can't be calculated at +/- epsilon,
      # happens when theta is at the edge of allowing matrix to be Cholesky
      # decomposed. Check here for that, use one side approx if only one grad
      # can be calculated. If both fail, no clue what to do.
      g_plus <- (CGGP_internal_gneglogpost(thetav,CGGP,y.thisloop,
                                           Xs=Xs, ys=ys.thisloop,
                                           HandlingSuppData=HandlingSuppData) *
                   (2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)/2
      
      g_minus <- (CGGP_internal_gneglogpost(thetav2,CGGP,y.thisloop,
                                            Xs=Xs, ys=ys.thisloop,
                                            HandlingSuppData=HandlingSuppData) *
                    (2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)/2
      
      if (all(is.finite(g_plus)) && all(is.finite(g_minus))) {
        H[ipara,] <- g_plus - g_minus
      } else {
        dimensions_that_need_fixing <- c(dimensions_that_need_fixing, ipara)
        # message(c("At least one was not finite, ", g_plus, g_minus))
        if (all(is.finite(g_plus))) {
          H[ipara,] <- 2 * g_plus
        } else if (all(is.finite(g_minus))) {
          H[ipara,] <- -2 * g_minus
        } else {
          # stop("Having to set one to NaN, will probably break stuff")
          H[ipara,] <- NaN
        }
      }
      
    }
    
    # For any dimensions that gave issues, set them to unit vector here
    # Shouldn't affect eigen stuff for other dimensions?
    for (ipara in dimensions_that_need_fixing) {
      message(paste(c("Had to fix dimensions ", dimensions_that_need_fixing, " in CGGP_fit")))
      H[ipara,] <- H[,ipara] <- 0
      H[ipara, ipara] <- 1
    }
    
    Hmat = H/2+t(H)/2
    A = eigen(Hmat)
    
    cHa = (A$vectors)%*%diag(pmin(sqrt(pmax(1/(A$values),10^(-16))),4))%*%t(A$vectors)
    
    # Get posterior samples using Laplace approximation
    PST= log((1+thetaMAP)/(1-thetaMAP)) +
      cHa%*%matrix(rnorm(CGGP$numPostSamples*length(thetaMAP),0,1),
                   nrow=length(thetaMAP))
    thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    
    # Now if there were any bad dimensions, we need to set all those
    # thetaPostSamples for that dimension to be the MAP
    for (ipara in dimensions_that_need_fixing) {
      message(paste(c("Changed thetaPostSamples for dims ", dimensions_that_need_fixing)))
      thetaPostSamples[ipara,] <- thetaMAP[ipara]
    }
    
    if(CGGP$supplemented){
      # Cs = matrix(1,dim(CGGP$Xs)[1],CGGP$ss)
      # for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      #   V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$xb,
      #                    thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
      #   Cs = Cs*V[,CGGP$designindex[,dimlcv]]
      # }
      # 
      # Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
      # for (dimlcv in 1:CGGP$d) { # Loop over dimensions
      #   V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv],
      #                    thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
      #   Sigma_t = Sigma_t*V
      # }
      # 
      # MSE_s = list(matrix(0,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1]),
      #              (CGGP$d+1)*(CGGP$maxlevel+1)) 
      # for (dimlcv in 1:CGGP$d) {
      #   for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      #     MSE_s[[(dimlcv)*CGGP$maxlevel+levellcv]] =
      #       (-CGGP_internal_postvarmatcalc(CGGP$Xs[,dimlcv],CGGP$Xs[,dimlcv],
      #                                      CGGP$xb[1:CGGP$sizest[levellcv]],
      #                                      thetaMAP[(dimlcv-1)*CGGP$numpara +
      #                                                 1:CGGP$numpara],
      #                                      CorrMat=CGGP$CorrMat))
      #   }
      # }
      # 
      # for (blocklcv in 1:CGGP$uoCOUNT) {
      #   ME_s = matrix(1,nrow=dim(Xs)[1],ncol=dim(Xs)[1])
      #   for (dimlcv in 1:CGGP$d) {
      #     levelnow = CGGP$uo[blocklcv,dimlcv]
      #     ME_s = ME_s*MSE_s[[(dimlcv)*CGGP$maxlevel+levelnow]]
      #   }
      #   Sigma_t = Sigma_t-CGGP$w[blocklcv]*(ME_s)
      # }
      # yhats = Cs%*%pw
      # 
      # 
      # Sti_resid = solve(Sigma_t,ys.thisloop-yhats)
      # Sti = solve(Sigma_t)
      # sigma2MAP = (sigma2MAP*dim(CGGP$design)[1] +
      #                colSums((ys.thisloop-yhats)*Sti_resid)) / (
      #                  dim(CGGP$design)[1]+dim(Xs)[1])
      # 
      # pw_adj_y = t(Cs)%*%Sti_resid
      # pw_adj <- CGGP_internal_calcpw(CGGP=CGGP, y=pw_adj_y, theta=thetaMAP)
      # 
      # pw_uppadj = pw-pw_adj
      # supppw = Sti_resid
      
      supp_values <- CGGP_internal_calc_supp_pw_sigma2_Sti(
        CGGP, thetaMAP=thetaMAP, ys.thisloop=ys.thisloop, pw=pw,
        sigma2MAP=sigma2MAP, only_sigma2MAP=FALSE)
      supppw <- supp_values$supppw
      sigma2MAP <- supp_values$sigma2MAP
      Sti <- supp_values$Sti
      pw_uppadj<- supp_values$pw_uppadj
    }
    
    
    # Add all new variables to CGGP that are needed
    if (nopd==1) { # Only 1 output parameter dim, so just set them
      CGGP$thetaMAP <- thetaMAP
      CGGP$sigma2MAP <- sigma2MAP
      CGGP$pw <- pw
      CGGP$thetaPostSamples <- thetaPostSamples
      CGGP$cholSs = cholSs
      if (CGGP$supplemented) {
        CGGP$pw_uppadj <- pw_uppadj
        CGGP$supppw <- supppw
        CGGP$Sti = Sti
        CGGP$sigma2MAP <- sigma2MAP
      }
    } else { # More than 1 opd, so need to set as columns of matrix
      if (opdlcv==1) { # First time, initialize matrix/array for all
        
        CGGP$thetaMAP <- matrix(NaN, length(thetaMAP), nopd)
        if (length(sigma2MAP) != 1) {
          stop("Error: sigma2map should be a 1x1 matrix.")
        }
        CGGP$sigma2MAP <- numeric(nopd)
        CGGP$pw <- matrix(NaN, length(pw), nopd) 
        # thetaPostSamples is matrix, so this is 3dim array below
        CGGP$thetaPostSamples <- array(data = NaN,
                                       dim=c(dim(thetaPostSamples), nopd))
        CGGP$cholSs <- vector("list", nopd)
      }
      CGGP$thetaMAP[,opdlcv] <- thetaMAP
      CGGP$sigma2MAP[opdlcv] <- sigma2MAP
      
      CGGP$pw[,opdlcv] <- pw
      CGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      # CGGP$cholSs[,,opdlcv] <- cholSs
      CGGP$cholSs[[opdlcv]] <- cholSs
      if (CGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          CGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nopd)
          CGGP$supppw <- matrix(NaN, nrow(supppw), nopd)
          # Sti is matrix, so this is 3 dim array
          CGGP$Sti = array(NaN, dim=c(dim(Sti), nopd))
        }
        CGGP$pw_uppadj[,opdlcv] <- pw_uppadj
        CGGP$supppw[,opdlcv] <- supppw
        CGGP$Sti[,,opdlcv] = Sti
      }
    }
  }
  
  # Clear old sigma2_samples. They will be recalculated when needed.
  CGGP$sigma2_samples <- NULL
  
  return(CGGP)
}



#' Calculate posterior variance
#'
#' @param x1 Points at which to calculate MSE
#' @param x2 Levels along dimension, vector???
#' @param xo No clue what this is
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix
#' for vector of 1D points.
#' @param returndPVMC Should dPVMC be returned?
#' @param returndiagonly Should only the diagonal be returned?
#'
#' @return Variance posterior
# @export
#' @noRd
#'
#' @examples
#' CGGP_internal_postvarmatcalc(c(.4,.52), c(0,.25,.5,.75,1),
#'              xo=c(.11), theta=c(.1,.2,.3),
#'              CorrMat=CGGP_internal_CorrMatCauchySQT)
CGGP_internal_postvarmatcalc <- function(x1, x2, xo, theta, CorrMat,
                                         returndPVMC = FALSE,
                                         returndiagonly=FALSE) {
  if(!returndiagonly && !returndPVMC){
    S = CorrMat(xo, xo, theta)
    n = length(xo)
    cholS = chol(S)
    
    C1o = CorrMat(x1, xo, theta)
    CoinvC1o = backsolve(cholS,backsolve(cholS,t(C1o), transpose = TRUE))
    C2o = CorrMat(x2, xo, theta)
    Sigma_mat = - t(CoinvC1o)%*%t(C2o)  
    return(Sigma_mat)
  } else {
    stop("Full postvarmatcalc function was removed #25082")
    # Only the chunk above was ever used in our code,
    # the full version where the other options can be used
    # was moved to scratch/scratch_postvarmatcalc_fullversion.R
  }
}


#' Calculate sigma2 for all theta samples
#'
#' @param CGGP CGGP object
#'
#' @return All sigma2 samples
## @export
#' @noRd
CGGP_internal_calc_sigma2_samples <- function(CGGP) {
  nopd <- if (is.matrix(CGGP$thetaPostSamples)) {1} else {dim(CGGP$thetaPostSamples)[3]}
  
  if (is.null(CGGP[["y"]]) || length(CGGP$y)==0) { # Only supp data
    if (nopd == 1 && length(CGGP$sigma2MAP)==1) { # 1 opd and 1 od
      # Single output dimension
      as.matrix(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calc_supp_only_supppw_sigma2_Sti(
                  CGGP=CGGP,thetaMAP=th,ys.thisloop=CGGP$ys, only_sigma2MAP=TRUE
                )$sigma2
              }
        )
      )
    } else if (nopd == 1) { # 1 opd but 2+ od
      # MV output but shared parameters, so sigma2 is vector
      t(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calc_supp_only_supppw_sigma2_Sti(
                  CGGP=CGGP,thetaMAP=th,ys.thisloop=CGGP$ys, only_sigma2MAP=TRUE
                )$sigma2
              }
        )
      )
    } else { # 2+ opd, so must be 2+ od
      # MV output with separate parameters, so need to loop over
      #  both samples and output dimension
      outer(1:CGGP$numPostSamples, 1:nopd,
            Vectorize(
              function(samplenum, outputdim) {
                CGGP_internal_calc_supp_only_supppw_sigma2_Sti(
                  CGGP=CGGP,thetaMAP=CGGP$thetaPostSamples[,samplenum,outputdim],
                  ys.thisloop=CGGP$ys[,outputdim],
                  only_sigma2MAP=TRUE
                )$sigma2
              }
            )
      )
    }
    
    
  } else if (!CGGP$supplemented) {
    if (nopd == 1 && length(CGGP$sigma2MAP)==1) { # 1 opd and 1 od
      # Single output dimension
      as.matrix(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calcsigma2(CGGP,
                                         CGGP$y,
                                         th
                )$sigma2
              }
        )
      )
    } else if (nopd == 1) { # 1 opd but 2+ od
      # MV output but shared parameters, so sigma2 is vector
      t(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calcsigma2(CGGP,
                                         CGGP$y,
                                         th
                )$sigma2
              }
        )
      )
    } else { # 2+ opd, so must be 2+ od
      # MV output with separate parameters, so need to loop over
      #  both samples and output dimension
      outer(1:CGGP$numPostSamples, 1:nopd,
            Vectorize(function(samplenum, outputdim) {
              CGGP_internal_calcsigma2(
                CGGP,
                CGGP$y[,outputdim],
                CGGP$thetaPostSamples[,samplenum,outputdim]
              )$sigma2
            })
      )
    }
  } else { # There is supplementary data
    if (nopd == 1 && length(CGGP$sigma2MAP)==1) { # 1 opd and 1 od
      # Single output dimension
      as.matrix(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calc_supp_pw_sigma2_Sti(CGGP=CGGP,
                                                      y.thisloop=CGGP$y,
                                                      thetaMAP=th,
                                                      ys.thisloop=CGGP$ys,
                                                      only_sigma2MAP=TRUE
                )$sigma2
              }
        )
      )
    } else if (nopd == 1) { # 1 opd but 2+ od
      # MV output but shared parameters, so sigma2 is vector
      t(
        apply(CGGP$thetaPostSamples, 2,
              function(th) {
                CGGP_internal_calc_supp_pw_sigma2_Sti(CGGP,
                                                      y.thisloop=CGGP$y,
                                                      thetaMAP=th,
                                                      ys.thisloop=CGGP$ys,
                                                      only_sigma2MAP=TRUE
                )$sigma2
              }
        )
      )
    } else { # 2+ opd, so must be 2+ od
      # MV output with separate parameters, so need to loop over
      #  both samples and output dimension
      outer(1:CGGP$numPostSamples, 1:nopd,
            Vectorize(function(samplenum, outputdim) {
              CGGP_internal_calc_supp_pw_sigma2_Sti(
                CGGP=CGGP,
                y.thisloop=CGGP$y[,outputdim],
                thetaMAP=CGGP$thetaPostSamples[,samplenum,outputdim],
                ys.thisloop=CGGP$ys[,outputdim],
                only_sigma2MAP=TRUE
              )$sigma2
            })
      )
    }
  }
}
