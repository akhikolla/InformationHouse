#' Raman spectrum of methanol (CH3OH)
#' 
#' @format A list containing 2 variables:
#' \describe{
#'   \item{wavenumbers}{a numeric Vector of 331 wavenumbers (cm^-1)}
#'   \item{wavenumbers}{a \code{1 * 331} Matrix of intensity values (a.u.)}
#'  }
"methanol"

#' Surface-enhanced Raman spectram of tetramethylrhodamine+DNA (T20)
#' 
#' @format A list containing 2 variables:
#' \describe{
#'   \item{wavenumbers}{a numeric Vector of 2401 wavenumbers (cm^-1)}
#'   \item{wavenumbers}{a \code{1 * 2401} Matrix of intensity values (a.u.)}
#'  }
"lsTamra"

#' SMC particles for TAMRA+DNA (T20)
#'
#' Posterior distribution for pseudo-Voigt parameters, obtained by running
#' `fitVoigtPeaksSMC` on a spectrum from Gracie et al. (Anal. Chem., 2016).
#' 1000 SMC particles with 32 peaks. For details, see the vignette.
#'
#' @format A list containing 15 variables:
#' \describe{
#'   \item{weights}{normalised importance weights for each particle}
#'   \item{location}{location parameters of 32 peaks}
#'   \item{beta}{amplitudes of 32 peaks}
#'   \item{scale_G}{scale of the Gaussian (RBF) broadening}
#'   \item{scale_L}{scale of the Lorentzian (Cauchy) broadening}
#'   \item{sigma}{standard deviation of the additive white noise}
#'   \item{lambda}{smoothing parameter of the cubic B-splines}
#'   \item{priors}{List of informative priors}
#'   \item{ess}{history of the effective sample size}
#'   \item{kappa}{history of the likelihood tempering}
#'   \item{accept}{history of Metropolis-Hastings acceptance rates}
#'   \item{mhSteps}{history of Metropolis-Hastings steps}
#'   \item{times}{history of times for each SMC iteration}
#'   \item{time}{computation time taken by the SMC algorithm}
#' }
"result"

#' SMC particles for methanol (CH3OH)
#'
#' Posterior distribution for pseudo-Voigt parameters, obtained by running
#' `fitVoigtPeaksSMC` on a Raman spectrum of methanol with 4 peaks.
#' For details, refer to the vignette.
#'
#' @format A list containing 15 variables.
"result2"
