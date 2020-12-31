#' Miscellaneous functions to support the \code{ibmcraftr} packare are here.
#'
#' @param rates A numeric scalar or vector to be transformed into rates.
#' @return A numeric scalar or vector in terms of probabilities.
#' @examples
#' rate2prob(c(.1, .5))
#'
#' @export
#'

#misc functions

rate2prob <- function(rates){
  1-exp(-rates*1)
}

#unload
.onUnload <- function(libpath) { library.dynam.unload("ibmcraftr", libpath) }
