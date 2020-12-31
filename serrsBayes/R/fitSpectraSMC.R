#' Fit the model using Sequential Monte Carlo (SMC).
#'
#' @param wl Vector of \code{nwl} wavenumbers at which the spetra are observed.
#' @param spc \code{n_y * nwl} Matrix of observed Raman spectra.
#' @param peakWL Vector of locations for each peak (cm^-1)
#' @param lPriors List of hyperparameters for the prior distributions.
#' @param conc Vector of \code{n_y} nanomolar (nM) dye concentrations for each observation.
#' @param npart number of SMC particles to use for the importance sampling distribution.
#' @param rate the target rate of reduction in the effective sample size (ESS).
#' @param minESS minimum effective sample size, below which the particles are resampled.
#' @param destDir destination directory to save intermediate results (for long-running computations)
#'
#' @return a List containing weighted parameter values, known as particles:
#' \describe{
#'   \item{\code{weights}}{Vector of importance weights for each particle.}
#'   \item{\code{beta}}{\code{npart * npeaks} Matrix of regression coefficients for the amplitudes.}
#'   \item{\code{scale}}{\code{npart * npeaks} Matrix of scale parameters.}
#'   \item{\code{sigma}}{Vector of \code{npart} standard deviations.}
#'   \item{\code{alpha}}{\code{bl.knots * n_y * npart} Array of spline coefficients for the baseline.}
#'   \item{\code{basis}}{A dense \code{nwl * bl.knots} Matrix containing the values of the basis functions.}
#'   \item{\code{expFn}}{\code{npart * nwl} Matrix containing the spectral signature.}
#'   \item{\code{ess}}{Vector containing the effective sample size (ESS) at each SMC iteration.}
#'   \item{\code{logEvidence}}{Vector containing the logarithm of the model evidence (marginal likelihood).}
#'   \item{\code{accept}}{Vector containing the Metropolis-Hastings acceptance rate at each SMC iteration.}
#'   \item{\code{sd.mh}}{\code{niter * 2npeaks} Matrix of random walk MH bandwidths at each SMC iteration..}
#'   }
#' @references
#' Chopin (2002) "A Sequential Particle Filter Method for Static Models," Biometrika 89(3): 539--551,
#' DOI: \href{http://dx.doi.org/10.1093/biomet/89.3.539}{10.1093/biomet/89.3.539}
#' @importFrom methods as
#' @importFrom stats rlnorm rnorm rgamma
#' @examples 
#' wavenumbers <- seq(200,600,by=10)
#' spectra <- matrix(nrow=1, ncol=length(wavenumbers))
#' peakLocations <- c(300,500)
#' peakAmplitude <- c(10000,4000)
#' peakScale <- c(10, 15)
#' signature <- weightedLorentzian(peakLocations, peakScale, peakAmplitude, wavenumbers)
#' baseline <- 1000*cos(wavenumbers/200) + 2*wavenumbers
#' spectra[1,] <- signature + baseline + rnorm(length(wavenumbers),0,200)
#' lPriors <- list(scale.mu=log(11.6) - (0.4^2)/2, scale.sd=0.4, bl.smooth=10^11, bl.knots=20,
#'                 beta.mu=5000, beta.sd=5000, noise.sd=200, noise.nu=4)
#' result <- fitSpectraSMC(wavenumbers, spectra, peakLocations, lPriors, npart=500)
fitSpectraSMC <- function(wl, spc, peakWL, lPriors, conc=rep(1.0,nrow(spc)), npart=10000, rate=0.9, minESS=npart/2, destDir=NA) {
  NP <- length(peakWL)
  nWL <- length(wl)
  niter <- y.n <- nrow(spc)
  lPriors$noise.SS <- lPriors$noise.nu * lPriors$noise.sd^2
  n_acc <- ess <- evidence <- wlIDX <- rep(0,niter*nWL)
  sd.mh <- matrix(NA, nrow=niter*nWL, ncol=2*NP)
  print(paste("SMC with",y.n,"observations at",length(unique(conc)),"unique concentrations,",npart,"particles, and",nWL,"wavenumbers."))
  
  # Step 0: cubic B-spline basis (Green & Silverman 1994, Sect. 2.3.3)
  ptm <- proc.time()
  basisFn <- getBsplineBasis(wl, lPriors$bl.knots, lPriors$bl.smooth)
  lPriors$bl.knots <- ncol(basisFn$basis)
  lPriors$bl.basis <- basisFn$basis
  lPriors$bl.precision <- as(basisFn$precision, "dgCMatrix") # cast to the parent class
  gi <- t(lPriors$bl.basis)%*%lPriors$bl.basis + lPriors$bl.precision
  Sgi<-solve(gi,sparse=TRUE)
  bl.smoother <- tcrossprod(Sgi, lPriors$bl.basis)
  print(paste("Step 0: computing",lPriors$bl.knots,"B-spline basis functions took",(proc.time() - ptm)[3],"sec."))
  
  # Step 1: Initialization (draw particles from the prior)
  ptm <- proc.time()
  theta <- list()
  theta$beta <- theta$scale <- matrix(NA, ncol=npart, nrow=NP)
  theta$expFn <- matrix(0, ncol=npart, nrow=nWL)
  theta$alpha <- array(NA, dim=c(lPriors$bl.knots, y.n, npart))
  theta$sigma <- 1/sqrt(rgamma(npart, shape=lPriors$noise.nu/2, rate=lPriors$noise.SS/2))
  for (pt in 1:npart) {
    theta$scale[,pt] <- rlnorm(NP, meanlog=lPriors$scale.mu, sdlog=lPriors$scale.sd)
    ampMx <- rnorm(NP, lPriors$beta.mu, lPriors$beta.sd)
    ampMx[ampMx<0] <- lPriors$beta.mu # enforce the positivity constraint
    theta$beta[,pt] <- ampMx
    if ("peaks" %in% names(lPriors) && !is.null(lPriors$peaks) && lPriors$peaks != "Lorentzian") {
      theta$expFn[,pt] <- weightedGaussian(peakWL, theta$scale[,pt], theta$beta[,pt], wl)
    } else {
      theta$expFn[,pt] <- weightedLorentzian(peakWL, theta$scale[,pt], theta$beta[,pt], wl)
    }
  }
  logWt <- rep(log(1/npart), npart)
  ess[1] <- npart
  print(paste("Step 1: initialization for",NP,"peaks took",(proc.time() - ptm)[3],"sec."))

  it <- 0
  for (i in 1:y.n) {
    print(paste("Observation", i, format(Sys.time(), "%Y%m%d-%H%M%S")))
    ptm <- proc.time()
    j <- 1
    # initialise baseline for the new observation
    for (pt in 1:npart) {
      Obsi <- spc[i,] - conc[i]*theta$expFn[,pt]
      theta$alpha[,i,pt] <- as.vector(bl.smoother %*% Obsi)
    }
    print(paste("Baselines took",(proc.time() - ptm)[3],"sec"))
    reWtIDX <- sample.int(nWL)
    while (j <= length(wl)) {
      it <- it+1

      ptm1 <- proc.time()
      # Step 2: Adaptation (calculate importance weights)
      bl.est <- lPriors$bl.basis %*% theta$alpha[,i,]
      newWeights <- reWeightParticles(spc, conc[i]*theta$expFn, bl.est, i, j, theta$sigma, logWt, rate, reWtIDX)
      logWt <- newWeights$weights
      ess[it] <- newWeights$ess
      j <- newWeights$index+1
      evidence[it] <- newWeights$evidence
      # importance sampling estimate of the variance (for RWMH proposal bandwidth)
      beta.mu <- weightedMean(theta$beta, logWt)
      beta.sd <- 0.5*sqrt(2*weightedVariance(theta$beta, logWt, beta.mu)/sqrt(2*NP))
      scale.mu <- weightedMean(theta$scale, logWt)
      scale.sd <- 0.5*sqrt(2*weightedVariance(theta$scale, logWt, scale.mu)/sqrt(2*NP))
      sd.mh[it,] <- c(beta.sd, scale.sd)
      ptm4 <- proc.time()
      print(paste("reweighting took",(ptm4 - ptm1)[3],"sec for ESS",ess[it]))
      
      # Step 3: Resampling
      if (ess[it] < minESS && j < length(wl)) {
        idx <- resampleParticles(logWt, theta$beta, theta$scale, theta$expFn, theta$alpha, y.n, lPriors$bl.knots)
        theta$sigma <- theta$sigma[idx]
        logWt <- rep(log(1/npart), npart)
        ess[it] <- npart
        print(paste("*** Resampling with",length(unique(idx)),"unique indices took", (proc.time() - ptm4)[3],"sec ***"))
      }
      
      # Step 4: Mutation (Metropolis-Hastings step)
      ptm2 <- proc.time()
      n_acc[it] <- marginalMetropolisUpdate(spc, i, conc[1:i], wl, peakWL, theta$beta, theta$scale,
                                            theta$sigma, theta$expFn, theta$alpha, sd.mh[it,], lPriors)
      print(paste("MH acceptance rate",n_acc[it]/npart,"took",(proc.time() - ptm2)[3],"sec."))
      if (n_acc[it]/npart < 0.01) break;
    }
    gc()
    # save partial results to disk after every iteration
    theta$weights=exp(logWt)
    theta$logEvidence=evidence[1:it]
    theta$accept=n_acc[1:it]/npart
    theta$ess=ess[1:it]
    theta$sd.mh=sd.mh[1:it,]
    if (!is.na(destDir) && file.exists(destDir)) {
      fname <- file.path(dirname(destDir),basename(destDir),paste0("smcRaman",i,"of",y.n,"W",nWL,"Q",npart,".rda"))
      save(theta, file=fname)
      print(paste("Interim results saved to",fname))
    }
    print(paste("SMC iter.",i,"took",(proc.time() - ptm)[3],"sec."))
  }
  theta$beta <- t(theta$beta)
  theta$scale <- t(theta$scale)
  theta$expFn <- t(theta$expFn)
  theta$basis <- lPriors$bl.basis
  return(theta)
}
