#' Fit CGGP to only supplemental data
#' 
#' If only supplemental data is given, and there is no grid data,
#' then this function is used.
#' It basically does the standard GP calculations since there is
#' nothing from the grid.
#'
#' @param CGGP CGGP object
#' @param Xs X supp matrix
#' @param Ys Y supp data
#' @param separateoutputparameterdimensions Should output dimensions be fit separately?
#' @param theta0 Initial theta0 for optimization
#' @param numPostSamples How many posterior samples should be calculated?
#' @param set_thetaMAP_to Value for thetaMAP to be set to
#'
#' @return CGGP object
## @export
#' @noRd
#' @examples
#' d <- 3
#' n <- 30
#' Xs <- matrix(runif(d*n), n, d)
#' Ys <- apply(Xs, 1, function(x){x[1]/(x[3]+.2)+exp(x[3])*cos(x[2]^2)})
#' cg <- CGGPcreate(d, Xs=Xs, Ys=Ys, batchsize=0)
CGGP_internal_fitwithonlysupp <- function(CGGP, Xs, Ys,
                                          separateoutputparameterdimensions=TRUE,
                                          theta0=rep(0,CGGP$numpara*CGGP$d),
                                          set_thetaMAP_to,
                                          numPostSamples=NULL
) {
  # ==========================.
  # ====    Set Values    ====
  # ==========================.
  CGGP$supplemented = TRUE
  CGGP$Xs = Xs
  CGGP$Ys = Ys
  if (!is.null(numPostSamples)) {
    CGGP$numPostSamples <- numPostSamples
  }
  rm(numPostSamples)
  
  if(!is.matrix(Ys)){ # 1D output
    CGGP$mu = mean(Ys)
    # y = Y-CGGP$mu
    ys = Ys-CGGP$mu
  } else{ # Multiple outputs
    CGGP$mu = colMeans(Ys) # Could use Y, or colMeans(rbind(Y, Ys)), or make sure Ys is big enough for this
    # y <- sweep(Y, 2, CGGP$mu)
    ys <- sweep(Ys, 2, CGGP$mu)
  }
  # CGGP$y = y
  CGGP$ys = ys
  
  
  
  
  
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions && is.matrix(ys)) {
    ncol(ys)
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
    if (ncol(theta0) > nopd) {
      theta0 <- theta0[,1:nopd]
    } else {
      theta0 <- cbind(theta0, matrix(0,nrow(theta0), nopd-ncol(theta0)))
    }
  }
  
  
  # =================================================================.
  # ==== Optimize parameters for each output parameter dimension ====
  # =================================================================.
  
  for (opdlcv in 1:nopd) { # output parameter dimension
    
    ys.thisloop <- if (nopd==1) {ys} else {ys[,opdlcv]}
    theta0.thisloop <- if (nopd==1) {theta0} else {theta0[,opdlcv]}
    
    if (!missing(set_thetaMAP_to) && !is.null(set_thetaMAP_to)) {
      opt.out <- list(par = if (nopd>1) {set_thetaMAP_to[,opdlcv]} else {set_thetaMAP_to})
    } else {
      # Sometimes it starts in a bad spot
      # Catch it with try, and start it at zero.
      # If that fails, we don't know what to do.
      neglogpost_par <- CGGP_internal_neglogpost(theta = theta0.thisloop,
                                                 CGGP=CGGP, y=NULL, Xs=Xs,
                                       ys=ys.thisloop)
      
      if (is.infinite(neglogpost_par)) {
        theta0_touse <- rep(0, CGGP$d)
      } else {
        theta0_touse <- theta0.thisloop
      }
      opt.out <- nlminb(
        theta0_touse,
        objective = CGGP_internal_neglogpost,
        gradient = CGGP_internal_gneglogpost,
        y = NULL,
        HandlingSuppData="Only",
        Xs=Xs,
        ys=ys.thisloop,
        CGGP = CGGP,
        control = list(rel.tol = 1e-8,iter.max = 500)
      )
    }
    
    # Set new theta
    thetaMAP <- opt.out$par
    
    # =====================================.
    # ====    Get posterior samples    ====
    # =====================================.
    
    totnumpara = length(thetaMAP)
    
    # H is the Hessian at thetaMAP with reverse transformation
    H = matrix(0,nrow=totnumpara,ncol=totnumpara)
    # Transform so instead of -1 to 1 it is -Inf to Inf. Mostly in -5 to 5 though.
    PSTn=  log((1+thetaMAP)/(1-thetaMAP))
    # Reverse transformation
    thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
    # Grad of reverse transformation function
    grad0 = CGGP_internal_gneglogpost(thetav,CGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)
    for(c in 1:totnumpara){
      rsad = rep(0,totnumpara)
      rsad[c] =10^(-3)
      PSTn=  log((1+thetaMAP)/(1-thetaMAP)) + rsad
      thetav=(exp(PSTn)-1)/(exp(PSTn)+1)
      H[c,] = (CGGP_internal_gneglogpost(thetav,CGGP,y=NULL, Xs=Xs, ys=ys.thisloop, HandlingSuppData="Only")*(2*(exp(PSTn))/(exp(PSTn)+1)^2)-grad0 )*10^(3)
    }
    Hmat = H/2+t(H)/2
    A = eigen(Hmat)
    cHa = (A$vectors)%*%diag(abs(A$values)^(-1/2))%*%t(A$vectors)
    
    # Get posterior samples using Laplace approx
    PST= log((1+thetaMAP)/(1-thetaMAP)) + cHa%*%matrix(rnorm(CGGP$numPostSamples*length(thetaMAP),0,1),nrow=length(thetaMAP))
    thetaPostSamples = (exp(PST)-1)/(exp(PST)+1)
    
    # =================================.
    # ====    Calculate supp stuff ====
    # # =================================.
    # Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
    # for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    #   V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    #   Sigma_t = Sigma_t*V
    # }
    # 
    # Sti_chol <- chol(Sigma_t + diag(CGGP$nugget, nrow(Sigma_t), ncol(Sigma_t)))
    # Sti <- chol2inv(Sti_chol)
    # 
    # # Use backsolve for stability
    # # supppw <- Sti %*% ys.thisloop
    # # supppw is same as Sti_resid in other files
    # supppw <- backsolve(Sti_chol, backsolve(Sti_chol, ys.thisloop, transpose = T))
    # if (is.matrix(supppw) && ncol(supppw)==1) {supppw <- as.vector(supppw)}
    # 
    # # sigma2MAP <- (t(ys.thisloop) %*% supppw) / nrow(Xs)
    # # if (is.matrix(sigma2MAP)) {sigma2MAP <- diag(sigma2MAP)}
    # sigma2MAP = colSums(as.matrix(ys.thisloop)*as.matrix(supppw))/nrow(Xs)
    supp_quantities <- CGGP_internal_calc_supp_only_supppw_sigma2_Sti(
      CGGP=CGGP,thetaMAP=thetaMAP,ys.thisloop=ys.thisloop, only_sigma2MAP=FALSE
    )
    supppw <- supp_quantities$supppw
    sigma2MAP <- supp_quantities$sigma2MAP
    Sti <- supp_quantities$Sti
    
    
    # Add all new variables to CGGP that are needed
    if (nopd==1) { # Only 1 output parameter dim, so just set them
      CGGP$thetaMAP <- thetaMAP
      CGGP$sigma2MAP <- sigma2MAP
      # CGGP$pw <- pw
      CGGP$thetaPostSamples <- thetaPostSamples
      if (CGGP$supplemented) {
        # CGGP$pw_uppadj <- pw_uppadj
        CGGP$supppw <- supppw
        CGGP$Sti = Sti
        CGGP$sigma2MAP <- sigma2MAP
      }
    } else { # More than 1 opd, so need to set as columns of matrix
      if (opdlcv==1) { # First time, initialize matrix/array for all
        
        CGGP$thetaMAP <- matrix(NaN, length(thetaMAP), nopd)
        if (length(sigma2MAP) != 1) {
          stop("ERROR HERE, sigma2map can be matrix??? It is always a 1x1 matrix from what I've seen before.")
        }
        CGGP$sigma2MAP <- numeric(nopd)
        # CGGP$pw <- matrix(NaN, length(pw), nopd) 
        # thetaPostSamples is matrix, so this is 3dim array below
        CGGP$thetaPostSamples <- array(data = NaN, dim=c(dim(thetaPostSamples), nopd))
      }
      CGGP$thetaMAP[,opdlcv] <- thetaMAP
      CGGP$sigma2MAP[opdlcv] <- sigma2MAP
      # CGGP$pw[,opdlcv] <- pw
      CGGP$thetaPostSamples[,,opdlcv] <- thetaPostSamples
      if (CGGP$supplemented) {
        if (opdlcv==1) { # First time initialize all
          
          # CGGP$pw_uppadj <- matrix(NaN, nrow(pw_uppadj), nopd)
          CGGP$supppw <- matrix(NaN, length(supppw), nopd) # Could be nrow supppw if 1 col matrix
          CGGP$Sti = array(NaN, dim=c(dim(Sti), nopd)) # Sti is matrix, so this is 3 dim array
        }
        # CGGP$pw_uppadj[,opdlcv] <- pw_uppadj
        CGGP$supppw[,opdlcv] <- supppw
        CGGP$Sti[,,opdlcv] = Sti
      }
    }
  }
  
  # Get sigma2 samples
  CGGP$sigma2_samples <- CGGP_internal_calc_sigma2_samples(CGGP)
  
  CGGP
}
