###
#' @title Rho (kappa)
#' 
#' @description 
#' This function calculates rho for an observed kappa value with associated set parameters
#' (testSetLength and OcSBaserate). Called by \code{\link{rho}}. A p-value is returned and if 
#' this value is less than 0.05, it is said that the handset does generalize to the entire set
#' 
#' @template rho-params
#' @param testSetLength The length of the \code{\link[=getTestSet]{testSet}} (ignored unless \emph{data} is an observed kappa value)
#' @param method set to "c" to calculate using the C++ implmentation. Defaults to "standard"
#' 
#' @export
#' @return rho for the given parameters
###
rhoK = function(
  x, 
  OcSBaserate,
  testSetLength,
  testSetBaserateInflation = 0, 
  OcSLength = 10000, 
  replicates = 800, 
  ScSKappaThreshold = 0.9, 
  ScSKappaMin = 0.40, 
  ScSPrecisionMin = 0.6, 
  ScSPrecisionMax = 1,
  method = "standard"
) {
  if (
    (x > 1 | x < 0) || 
    (OcSBaserate > 1 | OcSBaserate < 0)
  ) stop("x and OcSBaserate must be between 0 and 1.")

  if (OcSBaserate < 0 || ScSKappaThreshold < 0 || ScSKappaMin < 0) 
    stop("OcSBaserate, ScSKappaThreshold, and ScSKappaMin must be positive.")

  if (
    testSetLength %% 1 != 0 ||
    OcSLength %% 1 != 0 ||
    replicates %% 1 != 0
  ) stop("testSetLength, OcSLength, and replicates values must be a positive integer.")

  if (ScSKappaThreshold < ScSKappaMin) 
    stop("ScSKappaThreshold must be greater than ScSKappaMin.")

  if (
    (ScSPrecisionMin < 0 | ScSPrecisionMin > 1) ||
    (ScSPrecisionMax < 0 | ScSPrecisionMax > 1)
  ) stop("ScSPrecisionMin and ScSPrecisionMax must be between 0 and 1.")

  if (ScSPrecisionMax < ScSPrecisionMin) 
    stop("ScSPrecisionMax must be greater than ScSPrecisionMin.")

  result <- NULL
  if (method == "standard") {
    result <- calcRho(  x, OcSBaserate, testSetLength, testSetBaserateInflation, OcSLength, replicates, ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax)
  }
  else {
    result <- calcRho_c(x, OcSBaserate, testSetLength, testSetBaserateInflation, OcSLength, replicates, ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax)
  }
  return(result)
}
