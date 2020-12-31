###
#' @include rhoK.R
#' @title Rho (set)
#' 
#' @description 
#' This function calculates rho and kappa for a given \code{\link[=getTestSet]{testSet}}, and returns a list containing both values. Called by \code{\link{rho}}.
#' 
#' @export
#' 
#' @template rho-params
#' 
#' @return A list of the format:\describe{
#' \item{rho}{The rho of the \code{\link{codeSet}}}
#' \item{kappa}{The Cohen's Kappa of the \code{\link{codeSet}}}
#' }
###
rhoSet = function(
  x,
  OcSBaserate = NULL,
  testSetBaserateInflation = 0,
  OcSLength = 10000,
  replicates = 800,
  ScSKappaThreshold = 0.9,
  ScSKappaMin = 0.40,
  ScSPrecisionMin = 0.6,
  ScSPrecisionMax = 1
) {
  kappa = kappaSet(x);
  
  row.count = nrow(x);
  if(is.null(OcSBaserate)) {
    OcSBaserate = length(which(x[,1] == 1)) / row.count;
  }

  return(list(
    rho = rhoK(
      kappa, 
      OcSBaserate,
      row.count,
      testSetBaserateInflation,
      OcSLength,
      replicates,
      ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax
    ), 
    kappa = kappa,
    recall = rating_set_recall(x),
    precision = rating_set_precision(x)
  ))
}