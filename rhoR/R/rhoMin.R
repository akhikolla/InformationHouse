###
#' @title Rho Min
#' 
#' @description 
#' This function calculates the minimum testSetLength where it is possible to get a rho less than alpha for the given parameters of rho.
#' 
#' @param baserate A \code{\link{baserate}}
#' @param alpha The threshold of significance for rho (similar to an alpha level for a p value), defaulted to 0.05
#' @param inc An integer indicating by how much the testSetLength should increase each iteration
#' @param printInc A boolean indicating whether to print out each increment value with it's corresponding significance for rho
#' @param ... Any additional parameters passed into \code{\link{rho}}
#' 
#' @examples
#' #Add testSetBaserateInflation as an additional parameter
#' rhoMin(0.2, testSetBaserateInflation = 0.33)
#' 
#' #Add testSetBaserateInflation as well as changing inc and selecting printInc
#' rhoMin(0.2, inc = 5, printInc = TRUE, testSetBaserateInflation = 0.33)
#' 
#' @export
#' @return  The minimum length of testSet, to the nearest multiple of inc, greater than the minimum length, that would give a value where rho less than alpha becomes mathematically possible.
###
rhoMin = function(baserate, alpha = 0.05, inc = 10, printInc = FALSE, ...){
  rhoVal = 1
  testSetLength = 0
  if(inc%%1 != 0 | inc <= 0) stop("Inc value must be a positive integer.")
  if(alpha > 1 | alpha < 0) stop("Alpha must be between 0 and 1.")
  
  
  while(rhoVal > alpha){
    testSetLength = testSetLength + inc
    rhoVal = rho(1, baserate, testSetLength, ...)
    if(printInc) {
      cat(testSetLength, rhoVal,"\n")
    }
  }
  return(testSetLength)
  
}