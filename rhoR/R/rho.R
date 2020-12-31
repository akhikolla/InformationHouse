### 
#' @title Rho
#' 
#' @description 
#' This function calculates rho for a \code{\link[=getTestSet]{testSet}}, \code{\link{contingencyTable}}, or an observed kappa value with associated set parameters (testSetLength and OcSBaserate). 
#' 
#' @details 
#' Rho is a Monte Carlo rejective method of interrater reliability statistics, implemented here for Cohen's Kappa. Rho constructs a collection of data sets in which kappa is below a specified threshold, and computes the empirical distribution on kappa based on the specified sampling procedure. Rho returns the percent of the empirical distribution greater than or equal to an observed kappa. As a result, Rho quantifies the type 1 error in generalizing from an observed test set to a true value of agreement between two raters.
#' 
#' Rho starts with an observed kappa value, calculated on a subset of a \code{\link{codeSet}}, known as an observed \code{\link[=getTestSet]{testSet}}, and a \emph{kappa threshold} which indicates what is considered significant agreement between raters.
#' 
#' It then generates a collection of fully-coded, simulated 
#' \code{\link[=getTestSet]{codeSets}} (ScS), further described in 
#' \code{\link{createSimulatedCodeSet}}, all of which have a kappa value below the kappa 
#' threshold and similar properties as the original \code{\link{codeSet}}.
#' 
#' Then, kappa is calculated on a \code{\link[=getTestSet]{testSet}} sampled from each of the ScSs in the
#' collection to create a null hypothesis distribution. These \code{\link[=getTestSet]{testSets}} mirror the observed \code{\link[=getTestSet]{testSets}} in their size and sampling method.  How these \code{\link[=getTestSet]{testSets}} are sampled is futher described in \code{\link{getTestSet}}.
#' 
#' The null hypothesis is that the observed \code{\link[=getTestSet]{testSet}}, was sampled from a data set, which, if both raters were to code in its entirety, would result in a level of agreement below the kappa threshold.
#' 
#' For example, using an alpha level of 0.05, if the observed kappa is greater than 95 percent of the kappas in the null hypothesis distribution, the null hypothesis is rejected. Then one can conclude that the two raters would have acceptable agreement had they coded the entire data set.
#'
#' @template rho-params
#' @param testSetLength The length of the \code{\link[=getTestSet]{testSet}} (ignored unless \emph{data} is an observed kappa value)
#' 
#' @examples
#' # Given an observed kappa value
#' rho(x = 0.88, OcSBaserate = 0.2, testSetLength = 80)
#' 
#' # Given a test Set
#' rho(x = codeSet)
#' 
#' # Given a contingency Table
#' rho(x = contingencyTable)
#'
#' @export
#' @return rho and kappa for the given data and parameters (unless kappa is given)
###
rho = function(
  x,
  OcSBaserate = NULL,
  testSetLength = NULL,
  testSetBaserateInflation = 0,
  OcSLength = 10000,
  replicates = 800,
  ScSKappaThreshold = 0.9,
  ScSKappaMin = 0.40,
  ScSPrecisionMin = 0.6,
  ScSPrecisionMax = 1
) {

  ### Calculate rho given a kappa value
  if (is(x, "numeric")) {
    if (is.null(OcSBaserate)) { stop("Must give a baserate when a kappa value is given") }
    if (is.null(testSetLength)) { stop("Must give a testSetLength when a kappa value is given") }

    result <- rhoK(x, OcSBaserate,
               testSetLength,
               testSetBaserateInflation,
               OcSLength,
               replicates,
               ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax) 
  }

  ### Calculate rho given a code set
  else if (nrow(x) != 2) {
    result <- rhoSet(x, OcSBaserate,
                  testSetBaserateInflation,
                  OcSLength,
                  replicates,
                  ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax)
  }

  ### Calculate rho given a contingency table
  else {
    result <- rhoCT(x, OcSBaserate,
                 testSetBaserateInflation,
                 OcSLength,
                 replicates,
                 ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax)
  }

  return(result)
}