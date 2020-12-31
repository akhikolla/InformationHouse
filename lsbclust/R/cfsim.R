#' Compare Simulation Results
#' 
#' Generic function to compare simulation results in \pkg{lsbclust}.
#' 
#' @inheritParams cfsim.lsbclust
#' @seealso \code{\link{cfsim.lsbclust}}, \code{\link{cfsim.T3Clusf}}
#' @export
cfsim <- function(fitted, actual, method = c("diag", "cRand")) {
  UseMethod("cfsim", fitted)
}
