###
#' @title Rho using a file
#' 
#' @description 
#' This function calculates rho and kappa for a given \code{\link[=getTestSet]{testSet}} as defined by the file 
#' and columns (col1, col2), and returns a list containing both values. Called by \code{\link{rho}}.
#' 
#' @export
#' 
#' @param col1 The first column from file
#' @param col2 The second column from file
#' @template rho-params
#' 
#' @return A list of the format:\describe{
#' \item{rho}{The rho of the \code{\link{codeSet}}}
#' \item{kappa}{The Cohen's Kappa of the \code{\link{codeSet}}}
#' }
###
rho.file <- function(
  x, col1, col2, 
  OcSBaserate = NULL,
  testSetBaserateInflation = 0,
  OcSLength = 10000,
  replicates = 800,
  ScSKappaThreshold = 0.9,
  ScSKappaMin = .40,
  ScSPrecisionMin = 0.6,
  ScSPrecisionMax = 1
) {
  setFile = utils::read.csv(x)
  set = as.matrix(setFile[,c(col1, col2)])
  
  if (!any(is.na(set))) {
    rhoSet(
      set,
      OcSBaserate,
      testSetBaserateInflation, 
      OcSLength, replicates, 
      ScSKappaThreshold, ScSKappaMin, 
      ScSPrecisionMin, ScSPrecisionMax
    )
  }
  else {
    stop("The columns provided for col1 and col2 must not contain NA")
  }
}