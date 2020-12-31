###
#' @title Rho (contingency Table)
#' 
#' @description 
#' This function calculates rho and kappa for a given \code{\link{contingencyTable}}, and returns a list containing both values. Called by \code{\link{rho}}.
#' 
#' 
#' @template rho-params
#' 
#' @export
#' @return A list of the format:\describe{
#' \item{rho}{The rho of the \code{\link{contingencyTable}}}
#' \item{kappa}{The Cohen's Kappa of the \code{\link{contingencyTable}}}
#' }
###
rhoCT = function(
  x,
  OcSBaserate = NULL,
  testSetBaserateInflation = 0,
  OcSLength = 10000,
  replicates = 800,
  ScSKappaThreshold = 0.9,
  ScSKappaMin = .40,
  ScSPrecisionMin = 0.6,
  ScSPrecisionMax = 1
){
  if (any(x < 0)) {
    stop("Values in Contingency Table must be positive")
  }
  kappa = kappa_ct(x);
  
  if (is.null(OcSBaserate)) {
    OcSBaserate = (sum(x[1,])) / (sum(x))
  }
  
  rho_result <- rhoK(
    kappa, 
    OcSBaserate,
    sum(x),
    testSetBaserateInflation,
    OcSLength,
    replicates,
    ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax
  )
  
  return(list(
    rho = rho_result, 
    kappa = kappa,
    recall = rating_set_recall(x),
    precision = rating_set_precision(x)
  ))
}