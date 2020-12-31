#' Fit the model using Markov chain Monte Carlo.
#'
#' @param wl Vector of \code{nwl} wavenumbers at which the spetra are observed.
#' @param spc \code{n_y * nwl} Matrix of observed Raman spectra.
#' @param peakWL Vector of locations for each peak (cm^-1)
#' @param lPriors List of hyperparameters for the prior distributions.
#' @param sd_mh Vector of \code{2 * npeaks} bandwidths for the random walk proposals.
#' @param niter number of MCMC iterations per chain.
#' @param nchains number of concurrent MCMC chains.
#'
#' @return a List containing MCMC samples for the model parameters:
#' \describe{
#'   \item{\code{amplitude}}{\code{niter * nchains * npeaks} Array of amplitudes.}
#'   \item{\code{scale}}{\code{niter * nchains * npeaks} Array of scale parameters.}
#'   \item{\code{sigma}}{\code{niter * nchains} Matrix of standard deviations.}
#'   \item{\code{n_acc}}{The number of RWMH proposals that were accepted.}
#'   }
#' @seealso \code{\link{marginalMetropolisUpdate}}
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
#'                 amp.mu=5000, amp.sd=5000, noise.sd=200, noise.nu=4)
#' rw_bw <- c(100, 100, 2, 2)
#' result <- fitSpectraMCMC(wavenumbers, spectra, peakLocations, lPriors, rw_bw, 500)
#' result$n_acc
fitSpectraMCMC <- function(wl, spc, peakWL, lPriors, sd_mh, niter=10000, nchains=4) {
  y.n <- nrow(spc)
  NS <- length(wl)
  NP <- length(peakWL)
  samp.ampl <- samp.scale <- array(dim=c(niter,nchains,NP))
  samp.sigma <- matrix(nrow=niter,ncol=nchains)
  lPriors$noise.SS <- lPriors$noise.nu * lPriors$noise.sd^2
  lPriors$beta.mu <- lPriors$amp.mu
  lPriors$beta.sd <- lPriors$amp.sd
  n_acc <- 0
  print(paste("MCMC with",y.n,"observations of",NP,"peaks at",NS,"wavenumbers."))
  
  # Step 0: cubic B-spline basis (Green & Silverman 1994, Sect. 2.3.3)
  ptm <- proc.time()
  basisFn <- getBsplineBasis(wl, lPriors$bl.knots, lPriors$bl.smooth)
  lPriors$bl.knots <- ncol(basisFn$basis)
  lPriors$bl.basis <- basisFn$basis
  lPriors$bl.precision <- as(basisFn$precision, "dgCMatrix") # cast to the parent class
  gi <- t(lPriors$bl.basis)%*%lPriors$bl.basis + lPriors$bl.precision
  Sgi<-solve(gi,sparse=TRUE)
  print(paste("Step 0: computing",lPriors$bl.knots,"B-spline basis functions took",(proc.time() - ptm)[3],"sec."))
  
  # Step 1: Initialization (draw particles from the prior)
  ptm <- proc.time()
  theta <- list()
  theta$beta <- theta$scale <- matrix(NA, ncol=nchains, nrow=NP)
  theta$expFn <- matrix(0, ncol=nchains, nrow=NS)
  theta$alpha <- array(NA, dim=c(lPriors$bl.knots, y.n, nchains))
  theta$sigma <- vector(mode="numeric", length=nchains)
  for (pt in 1:nchains) {
    theta$scale[,pt] <- rlnorm(NP, meanlog=lPriors$scale.mu, sdlog=lPriors$scale.sd)
    ampMx <- rnorm(NP, lPriors$amp.mu, lPriors$amp.sd)
    ampMx[ampMx<0] <- lPriors$amp.mu
    theta$beta[,pt] <- ampMx
    if ("peaks" %in% names(lPriors) && !is.null(lPriors$peaks) && lPriors$peaks != "Lorentzian") {
      theta$expFn[,pt] <- weightedGaussian(peakWL, theta$scale[,pt], theta$beta[,pt], wl)
    } else {
      theta$expFn[,pt] <- weightedLorentzian(peakWL, theta$scale[,pt], theta$beta[,pt], wl)
    }
    for (i in 1:y.n) {
      Obsi<-spc[i,] - theta$expFn[,pt]
      mi<-Sgi%*%(t(lPriors$bl.basis)%*%Obsi)[,1]
      theta$alpha[,i,pt] <- as.vector(mi)
    }
  }
  theta$sigma <- 1/sqrt(rgamma(nchains, shape=lPriors$noise.nu/2, rate=lPriors$noise.SS/2))
  conc <- rep(1.0, y.n)
  print(paste("Step 1: initialization took",(proc.time() - ptm)[3],"sec for",nchains,"parallel chains."))
  
  n_acc <- 0
  for (it in 1:niter) {
    acc <- marginalMetropolisUpdate(spc, y.n, conc, wl, peakWL, theta$beta, theta$scale, theta$sigma, theta$expFn, theta$alpha, sd_mh, lPriors)
    samp.ampl[it,,] <- t(theta$beta)
    samp.scale[it,,] <- t(theta$scale)
    samp.sigma[it,] <- theta$sigma
    n_acc <- n_acc + acc
  }
  list(amplitude=samp.ampl, scale=samp.scale, sigma=samp.sigma, n_acc=n_acc, alpha=theta$alpha, peaks=theta$expFn, basis=lPriors$bl.basis)
}