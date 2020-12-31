###
#' @title Calculate kappa (Set)
#' 
#' @description 
#' This function calculates Cohen's kappa for a given \code{\link{codeSet}}. Called by \code{\link{kappa}}.
#' 
#' @export
#' @param set A \code{\link{codeSet}}
#' 
#' @seealso \code{\link{kappa}} and \code{\link{kappaCT}}
#' 
#' @return 
#' The kappa of the \code{\link{codeSet}}
#' 
###
kappaSet = function(set){
  if(any(set != 1 & set != 0)) {
    stop("Each set must only consist of zeros and ones.")
  }
  
  return(calcKappa(set,isSet = TRUE,kappaThreshold = NULL));
}

