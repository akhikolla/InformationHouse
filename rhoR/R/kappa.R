###
#' @title Calculate kappa
#'
#' @description
#' This function calculates Cohen's kappa on a \code{\link{contingencyTable}} or a \code{\link{codeSet}}
#'
#' @param data A \code{\link{contingencyTable}} or a \code{\link{codeSet}}
#'
#'
#' @seealso \code{\link{kappaSet}} and \code{\link{kappaCT}}
#'
#' @examples
#' #Given a code set
#' kappa(data = codeSet)
#'
#' #Given a contingency Table
#' kappa(data = contingencyTable)
#'
#' @export
#' @return The kappa of the \code{\link{contingencyTable}} or \code{\link{codeSet}}
###
kappa <- function(data) {
  if (all(dim(data) == 2)) {
    return(kappaCT(data))
  } else {
    return(kappaSet(data))
  }
}
