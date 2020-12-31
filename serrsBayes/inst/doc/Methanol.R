## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(1234)

## ----data, message=FALSE, fig.width=6-----------------------------------------
library(serrsBayes)
data("methanol", package = "serrsBayes")
wavenumbers <- methanol$wavenumbers
spectra <- methanol$spectra

peakLocations <- c(1033, 1106, 1149, 1448)
nPK <- length(peakLocations)
pkIdx <- numeric(nPK)
for (i in 1:nPK) {
  pkIdx[i] <- which.min(wavenumbers < peakLocations[i])
}
nWL <- length(wavenumbers)
plot(wavenumbers, spectra[1,], type='l', col=4,
     xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", main="Observed Raman spectrum for methanol")
points(peakLocations, spectra[1,pkIdx], cex=2, col=2)
text(peakLocations + c(100,20,40,0), spectra[1,pkIdx] + c(0,700,400,700), labels=peakLocations)

## ----fitSpectraSMC, results='hide'--------------------------------------------
lPriors <- list(scale.mu=log(11.6) - (0.4^2)/2, scale.sd=0.4, bl.smooth=10^11, bl.knots=50,
                 beta.mu=5000, beta.sd=5000, noise.sd=200, noise.nu=4)
tm <- system.time(result <- fitSpectraSMC(wavenumbers, spectra, peakLocations, lPriors))

## ----timeSpectra--------------------------------------------------------------
print(paste(tm["elapsed"]/60, "minutes"))

## ----fitSpectraBL, fig.width=6------------------------------------------------
samp.idx <- sample.int(length(result$weights), 200, prob=result$weights)
plot(wavenumbers, spectra[1,], type='l', lwd=3,
     xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", main="Fitted model with Lorentzian peaks")
for (pt in samp.idx) {
  bl.est <- result$basis %*% result$alpha[,1,pt]
  lines(wavenumbers, bl.est, col="#C3000009")
  lines(wavenumbers, bl.est + result$expFn[pt,], col="#00C30009")
}

## ----fitVoigtPeaksSMC, results='hide'-----------------------------------------
lPriors2 <- list(loc.mu=peakLocations, loc.sd=rep(50,nPK), scaG.mu=log(16.47) - (0.34^2)/2,
                 scaG.sd=0.34, scaL.mu=log(25.27) - (0.4^2)/2, scaL.sd=0.4, noise.nu=5,
                 noise.sd=50, bl.smooth=1, bl.knots=50)
data("result2", package = "serrsBayes")
if(!exists("result2")) {
  tm2 <- system.time(result2 <- fitVoigtPeaksSMC(wavenumbers, spectra, lPriors2, npart=3000))
  result2$time <- tm2
  save(result2, file="Figure 2/result2.rda")
}

## ----fitVoigtBL, fig.width=6--------------------------------------------------
print(paste(result2$time["elapsed"]/60, "minutes"))
samp.idx <- 1:nrow(result2$location)
plot(wavenumbers, spectra[1,], type='l', lwd=3,
     xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", main="Fitted model with Voigt peaks")
samp.mat <- resid.mat <- matrix(0,nrow=length(samp.idx), ncol=nWL)
samp.sigi <- samp.lambda <- numeric(length=nrow(samp.mat))
for (pt in 1:length(samp.idx)) {
  k <- samp.idx[pt]
  samp.mat[pt,] <- mixedVoigt(result2$location[k,], result2$scale_G[k,],
                              result2$scale_L[k,], result2$beta[k,], wavenumbers)
  samp.sigi[pt] <- result2$sigma[k]
  samp.lambda[pt] <- result2$lambda[k]
  
  Obsi <- spectra[1,] - samp.mat[pt,]
  g0_Cal <- length(Obsi) * samp.lambda[pt] * result2$priors$bl.precision
  gi_Cal <- crossprod(result2$priors$bl.basis) + g0_Cal
  mi_Cal <- as.vector(solve(gi_Cal, crossprod(result2$priors$bl.basis, Obsi)))

  bl.est <- result2$priors$bl.basis %*% mi_Cal # smoothed residuals = estimated basline
  lines(wavenumbers, bl.est, col="#C3000009")
  lines(wavenumbers, bl.est + samp.mat[pt,], col="#00C30009")
  resid.mat[pt,] <- Obsi - bl.est[,1]
}

## ----peakLoc, fig.cap="Posteriors for the peak locations.", fig.show='hold', fig.width=3----
for (j in 1:nPK) {
  hist(result2$location[,j], main=paste("Peak",j),
       xlab=expression(paste("Peak location (cm"^{-1}, ")")),
       freq=FALSE, xlim=range(result2$location[,j], peakLocations[j]))
  lines(density(result2$location[,j]), col=4, lty=2, lwd=3)
  abline(v=peakLocations[j], col=2, lty=3, lwd=3)
}

## ----voigt, fig.show='hold', fig.width=3--------------------------------------
nPart <- nrow(result2$beta)
result2$voigt <- result2$FWHM <- matrix(nrow=nPart, ncol=nPK)
for (k in 1:nPart) {
  result2$voigt[k,] <- getVoigtParam(result2$scale_G[k,], result2$scale_L[k,])
  f_G <- 2*result2$scale_G[k,]*sqrt(2*log(2))
  f_L <- 2*result2$scale_L[k,]
  result2$FWHM[k,] <- 0.5346*f_L + sqrt(0.2166*f_L^2 + f_G^2)
}
for (j in 1:nPK) {
  hist(result2$voigt[,j], main=paste("Peak",j),
       xlab=expression(eta), freq=FALSE, xlim=c(0,1))
  lines(density(result2$voigt[,j]), col=4, lty=2, lwd=3)
}

## ----fitSpectraSMCagain, results='hide'---------------------------------------
lPriors3 <- list(scale.mu=log(11.6) - (0.4^2)/2, scale.sd=0.4, bl.smooth=10^11, bl.knots=50, noise.sd=200, noise.nu=4)
lPriors3$beta.mu <- mean(result2$beta)
lPriors3$beta.sd <- sd(result2$beta)
pkLocNew <- c(colMeans(result2$location), 1550)
tm3 <- system.time(result3 <- fitSpectraSMC(wavenumbers, spectra, pkLocNew, lPriors3, npart=2000))

## ----time3--------------------------------------------------------------------
print(paste(tm3["elapsed"]/60, "minutes"))

## ----fitRedux, fig.width=6----------------------------------------------------
samp.idx <- sample.int(length(result3$weights), 200, prob=result3$weights)
plot(wavenumbers, spectra[1,], type='l', lwd=3,
     xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", main="Fitted model with informative priors")
for (pt in samp.idx) {
  bl.est <- result3$basis %*% result3$alpha[,1,pt]
  lines(wavenumbers, bl.est, col="#C3000009")
  lines(wavenumbers, bl.est + result3$expFn[pt,], col="#00C30009")
}
rug(pkLocNew, col=4, lwd=2)

