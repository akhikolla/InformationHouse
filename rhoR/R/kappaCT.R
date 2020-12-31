###
#' @title Calculate kappa (contingency Table)
#' 
#' @description 
#' This function calculates Cohen's kappa on a \code{\link{contingencyTable}}. Called by \code{\link{kappa}}. 
#' 
#' @param ct A \code{\link{contingencyTable}}
#' 
#' @seealso \code{\link{kappa}} and \code{\link{kappaSet}}
#' 
#' @export
#' @return  The kappa of the \code{\link{contingencyTable}}
###
kappaCT = function(ct){
  if(any(ct < 0)){stop("Values in Contingency Table must be positive")}
  if(ncol(ct) != 2 || nrow(ct) != 2) stop("Incorrect number of dimensions: Contingency table must be 2x2.")
  for(i in 1:nrow(ct)) {
    for(j in 1:ncol(ct)) {
      if (ct[i,j]%%1 != 0) stop("Contingency table values must be positive integers.")
    }  
  }
  
  return(calcKappa(ct,isSet = FALSE,kappaThreshold = NULL));
}
