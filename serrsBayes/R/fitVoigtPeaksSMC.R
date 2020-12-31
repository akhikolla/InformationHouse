#' Fit the model with Voigt peaks using Sequential Monte Carlo (SMC).
#'
#' @inheritParams fitSpectraSMC
#' @param mcAR target acceptance rate for the MCMC kernel
#' @param mcSteps number of iterations of the MCMC kernel
#' @importFrom methods as
#' @importFrom stats rlnorm rnorm rgamma runif cov.wt cov2cor median
#' @importFrom truncnorm rtruncnorm
#' @importFrom splines bs
#' @importFrom Matrix Matrix crossprod determinant
#' @examples 
#' wavenumbers <- seq(200,600,by=10)
#' spectra <- matrix(nrow=1, ncol=length(wavenumbers))
#' peakLocations <- c(300,500)
#' peakAmplitude <- c(10000,4000)
#' peakScale <- c(10, 15)
#' signature <- weightedLorentzian(peakLocations, peakScale, peakAmplitude, wavenumbers)
#' baseline <- 1000*cos(wavenumbers/200) + 2*wavenumbers
#' spectra[1,] <- signature + baseline + rnorm(length(wavenumbers),0,200)
#' lPriors <- list(scaG.mu=log(11.6) - (0.4^2)/2, scaG.sd=0.4, scaL.mu=log(11.6) - (0.4^2)/2,
#'    scaL.sd=0.4, bl.smooth=5, bl.knots=20, loc.mu=peakLocations, loc.sd=c(5,5),
#'    beta.mu=c(5000,5000), beta.sd=c(5000,5000), noise.sd=200, noise.nu=4)
#' result <- fitVoigtPeaksSMC(wavenumbers, spectra, lPriors, npart=50, mcSteps=1)
fitVoigtPeaksSMC <- function(wl, spc, lPriors, conc=rep(1.0,nrow(spc)), npart=10000, rate=0.9, mcAR=0.23, mcSteps=10, minESS=npart/2, destDir=NA) {
  N_Peaks <- length(lPriors$loc.mu)
  N_WN_Cal <- length(wl)
  N_Obs_Cal <- nrow(spc)
  lPriors$noise.SS <- lPriors$noise.nu * lPriors$noise.sd^2
  print(paste("SMC with",N_Obs_Cal,"observations at",length(unique(conc)),"unique concentrations,",npart,"particles, and",N_WN_Cal,"wavenumbers."))
  
  # Step 0: cubic B-spline basis (Eilers & Marx, 1996)
  ptm <- proc.time()
  Knots<-seq(min(wl),max(wl), length.out=lPriors$bl.knots)
  r <- max(diff(Knots))
  NK<-lPriors$bl.knots
  X_Cal <- bs(wl, knots=Knots, Boundary.knots = range(Knots) + c(-r,+r),
              intercept = TRUE)
  class(X_Cal) <- "matrix"
  XtX <- Matrix(crossprod(X_Cal), sparse=TRUE)
  NB_Cal<-ncol(X_Cal)
  FD2_Cal <- diff(diff(diag(NB_Cal))) # second derivatives of the spline coefficients
  Pre_Cal <- Matrix(crossprod(FD2_Cal), sparse=TRUE)
  
  # Demmler-Reinsch orthogonalisation
  R = chol(XtX + Pre_Cal*1e-9) # just in case XtX is ill-conditioned
  Rinv <- solve(R)
  Rsvd <- svd(crossprod(Rinv, Pre_Cal %*% Rinv))
  Ru <- Rinv %*% Rsvd$u
  A <- X_Cal %*% Rinv %*% Rsvd$u
  lPriors$bl.basis <- X_Cal
  lPriors$bl.precision <- as(Pre_Cal, "dgCMatrix")
  lPriors$bl.XtX <- as(XtX, "dgCMatrix")
  lPriors$bl.orthog <- as.matrix(A)
  lPriors$bl.Ru <- as.matrix(Ru)
  lPriors$bl.eigen <- Rsvd$d
  print(paste0("Step 0: computing ",NB_Cal," B-spline basis functions (r=",r,") took ",(proc.time() - ptm)[3],"sec."))

  # Step 1: Initialization (draw particles from the prior)
  ptm <- proc.time()
  Sample<-matrix(numeric(npart*(4*N_Peaks+3+N_Obs_Cal)),nrow=npart)
  Sample[,1:N_Peaks] <- rlnorm(N_Peaks*npart, lPriors$scaG.mu, lPriors$scaG.sd)
  Sample[,(N_Peaks+1):(2*N_Peaks)] <- rlnorm(N_Peaks*npart, lPriors$scaL.mu, lPriors$scaL.sd)
  # enforce window and identifiability constrants on peak locations
  for (k in 1:npart) {
    propLoc <- rtruncnorm(N_Peaks, a=min(wl), b=max(wl), mean=lPriors$loc.mu, sd=lPriors$loc.sd)
    Sample[k,(2*N_Peaks+1):(3*N_Peaks)] <- sort(propLoc)
  }
  # optional prior on beta
  if (exists("beta.mu", lPriors) && exists("beta.sd", lPriors)) {
    for (j in 1:N_Peaks) {
      Sample[,3*N_Peaks+j] <- rtruncnorm(npart, a=0, mean=lPriors$beta.mu[j], sd=lPriors$beta.sd[j])
    }
  } else { # otherwise, use uniform prior
    Sample[,(3*N_Peaks+1):(4*N_Peaks)] <- runif(N_Peaks*npart, 0, diff(range(spc))/max(conc))
  }
  Offset_1<-4*N_Peaks
  Offset_2<-Offset_1 + N_Obs_Cal + 1
  Cal_I <- 1
  Sample[,Offset_2+1] <- 1/rgamma(npart, lPriors$noise.nu/2, lPriors$noise.SS/2)
  #Sample[,Offset_1+4] <- 1/rgamma(npart, lPriors$bl.nu/2, lPriors$bl.SS/2)
  Sample[,Offset_2+2] <- Sample[,Offset_2+1]/lPriors$bl.smooth
  print(paste("Mean noise parameter sigma is now",mean(sqrt(Sample[,Offset_2+1]))))
  print(paste("Mean spline penalty lambda is now",mean(Sample[,Offset_2+1]/Sample[,Offset_2+2])))
  # compute log-likelihood:
  g0_Cal <- N_WN_Cal * lPriors$bl.smooth * Pre_Cal
  gi_Cal <- XtX + g0_Cal
  a0_Cal <- lPriors$noise.nu/2
  ai_Cal <- a0_Cal + N_WN_Cal/2
  b0_Cal <- lPriors$noise.SS/2
  for(k in 1:npart) {
        Sigi <- conc[Cal_I] * mixedVoigt(Sample[k,2*N_Peaks+(1:N_Peaks)], Sample[k,(1:N_Peaks)],
                                         Sample[k,N_Peaks+(1:N_Peaks)], Sample[k,3*N_Peaks+(1:N_Peaks)], wl)
        Obsi <- spc[Cal_I,] - Sigi
        lambda <- lPriors$bl.smooth # fixed smoothing penalty
        L_Ev <- computeLogLikelihood(Obsi, lambda, lPriors$noise.nu, lPriors$noise.SS,
                                          X_Cal, Rsvd$d, lPriors$bl.precision, lPriors$bl.XtX,
                                          lPriors$bl.orthog, lPriors$bl.Ru)
        Sample[k,Offset_1+2]<-L_Ev
  }
  Sample[,Offset_1+1]<-rep(1/npart,npart)
  T_Sample<-Sample
  T_Sample[,1:N_Peaks]<-log(T_Sample[,1:N_Peaks]) # scaG
  T_Sample[,(N_Peaks+1):(2*N_Peaks)]<-log(T_Sample[,(N_Peaks+1):(2*N_Peaks)]) # scaL
  T_Sample[,(3*N_Peaks+1):(4*N_Peaks)]<-log(T_Sample[,(3*N_Peaks+1):(4*N_Peaks)]) # amp/beta
  iTime <- proc.time() - ptm

  ESS<-1/sum(Sample[,Offset_1+1]^2)
  MC_Steps<-numeric(1000)
  MC_AR<-numeric(1000)
  ESS_Hist<-numeric(1000)
  ESS_AR<-numeric(1000)
  Kappa_Hist<-numeric(1000)
  Time_Hist<-numeric(1000)
  
  MC_Steps[1]<-0
  MC_AR[1]<-1
  ESS_Hist[1]<-ESS
  ESS_AR[1]<-npart
  Kappa_Hist[1]<-0
  Time_Hist[1]<-iTime[3]
  print(paste("Step 1: initialization for",N_Peaks,"Voigt peaks took",iTime[3],"sec."))
  print(colMeans(Sample[,(3*N_Peaks+1):(4*N_Peaks)]))

  i<-1
  Cal_I <- 1
  MADs<-numeric(4*N_Peaks)
  Alpha<-rate
  MC_AR[1]<-mcAR
  MCMC_MP<-1
  repeat{
    i<-i+1
    
    iTime<-system.time({

      ptm <- proc.time()
      Min_Kappa<-Kappa_Hist[i-1]
      Max_Kappa<-1
      Kappa<-1

      Temp_w<-Sample[,Offset_1+1]*exp((Kappa-Kappa_Hist[i-1])*(Sample[,Offset_1+Cal_I+1]-max(Sample[,Offset_1+Cal_I+1])))
      Temp_W<-Temp_w/sum(Temp_w)
      
      US1<-unique(Sample[,1])
      N_UP<-length(US1)
      
      Temp_W2<-numeric(N_UP)
      for(k in 1:N_UP){
        Temp_W2[k]<-sum(Temp_W[which(Sample[,1]==US1[k])])
      }
      
      Temp_ESS<-1/sum(Temp_W2^2)
      if(Temp_ESS<(Alpha*ESS_AR[i-1])){
        while(abs(Temp_ESS-((Alpha*ESS_AR[i-1])))>1 & !isTRUE(all.equal(Kappa, Min_Kappa))){
          if(Temp_ESS<((Alpha*ESS_AR[i-1]))){
            Max_Kappa<-Kappa
          } else{
            Min_Kappa<-Kappa
          }
          
          Kappa<-0.5*(Min_Kappa+Max_Kappa)
          
          Temp_w<-Sample[,Offset_1+1]*exp((Kappa-Kappa_Hist[i-1])*(Sample[,Offset_1+Cal_I+1]-max(Sample[,Offset_1+Cal_I+1])))
          Temp_W<-Temp_w/sum(Temp_w)
          
          US1<-unique(Sample[,1])
          N_UP<-length(US1)
          
          Temp_W2<-numeric(N_UP)
          for(k in 1:N_UP){
            Temp_W2[k]<-sum(Temp_W[which(Sample[,1]==US1[k])])
          }
          
          Temp_ESS<-1/sum(Temp_W2^2)
          
        }
      }
      
      Sample[,Offset_1+1]<-Temp_W
      Kappa_Hist[i]<-Kappa
      ESS_Hist[i]<-Temp_ESS
      
      print(paste0("Reweighting took ",(proc.time()-ptm)[3],"sec. for ESS ",Temp_ESS," with new kappa ",Kappa,"."))
      
      Acc<-0
      
      Prop_Info<-cov.wt(T_Sample[,1:(4*N_Peaks)],wt=Sample[,Offset_1+1])
      Prop_Mu<-Prop_Info$center
      Prop_Cor<-cov2cor(Prop_Info$cov)
      
      if(ESS_Hist[i] < minESS){
        # simple multinomial resampling
        ptm <- proc.time()
        ReSam<-sample(1:npart,size=npart,replace=T,prob=Sample[,Offset_1+1])
        Sample<-Sample[ReSam,]
        T_Sample<-T_Sample[ReSam,]
        
        Sample[,Offset_1+1]<-rep(1/npart,npart)
        T_Sample[,Offset_1+1]<-rep(1/npart,npart)
        print(paste("*** Resampling with",length(unique(T_Sample[,1])),"unique indices took",(proc.time()-ptm)[3],"sec ***"))
      }
      
      for(j in 1:(4*N_Peaks)){
        Prop_Mu[j]<-median(T_Sample[,j])
        MADs[j]<-median(abs((T_Sample[,j])-median(T_Sample[,j])))
      }
      
      Prop_Cov<-(1.4826*MADs)%*%t(1.4826*MADs)*Prop_Cor
      
      US1<-unique(T_Sample[,1])
      N_UP<-length(US1)
      
      Temp_W<-numeric(N_UP)
      for(k in 1:N_UP){
        Temp_W[k]<-sum(T_Sample[which(T_Sample[,1]==US1[k]),Offset_1+1])
      }
      
      Temp_ESS<-1/sum(Temp_W^2)
      ESS_AR[i]<-Temp_ESS
      
      if(!is.na(MC_AR[i-1])){
        MCMC_MP<-2^(-5*(0.23-MC_AR[i-1]))*MCMC_MP
      }
      mhCov <- MCMC_MP*(2.38^2/(4*N_Peaks))*Prop_Cov
      mhChol <- t(chol(mhCov, pivot = FALSE)) # error if not non-negative definite

      for(mcr in 1:mcSteps){
        MC_Steps[i]<-MC_Steps[i]+1
        mh_acc <- mhUpdateVoigt(spc, Cal_I, Kappa_Hist[i], conc, wl, Sample, T_Sample, mhChol, lPriors)
        Acc <- Acc + mh_acc

        # update effective sample size
        US1<-unique(Sample[,1])
        N_UP<-length(US1)
        Temp_W<-numeric(N_UP)
        for(k in 1:N_UP){
          Temp_W[k]<-sum(Sample[which(Sample[,1]==US1[k]),Offset_1+1])
        }
        Temp_ESS<-1/sum(Temp_W^2)
        print(paste(mh_acc,"M-H proposals accepted. Temp ESS is",Temp_ESS))
        ESS_AR[i]<-Temp_ESS
      }
      
      MC_AR[i]<-Acc/(npart*MC_Steps[i])
    })
    
    Time_Hist[i]<-iTime[3]
    
    if (!is.na(destDir) && file.exists(destDir)) {
      iFile<-paste0(destDir,"/Iteration_",i,"/")
      dir.create(iFile)
      save(Sample,file=paste0(iFile,"Sample.rda"))
      print(paste("Interim results saved to",iFile))
    }

    print(colMeans(Sample[,(3*N_Peaks+1):(4*N_Peaks)]))
    print(paste0("Iteration ",i," took ",iTime[3],"sec. for ",MC_Steps[i]," MCMC loops (acceptance rate ",MC_AR[i],")"))
    if (Kappa >= 1 || MC_AR[i] < 1/npart) {
      break
    }
  }

  if (Kappa < 1 && MC_AR[i] < 1/npart) {
    print(paste("SMC collapsed due to MH acceptance rate",
                Acc,"/",(npart*MC_Steps[i]),"=", MC_AR[i]))
  }
  return(list(priors=lPriors, ess=ESS_Hist[1:i], weights=Sample[,Offset_1+1], kappa=Kappa_Hist[1:i],
              accept=MC_AR[1:i], mhSteps=MC_Steps[1:i], essAR=ESS_AR[1:i], times=Time_Hist[1:i],
              scale_G=Sample[,1:N_Peaks], scale_L=Sample[,(N_Peaks+1):(2*N_Peaks)],
              location=Sample[,(2*N_Peaks+1):(3*N_Peaks)], beta=Sample[,(3*N_Peaks+1):(4*N_Peaks)],
              sigma=sqrt(Sample[,Offset_2+1]), lambda=Sample[,Offset_2+1]/Sample[,Offset_2+2]))
}
