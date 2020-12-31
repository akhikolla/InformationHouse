#' Print \code{evgram} object
#' 
#' Print an object of class \code{evgram} produced by the 
#' \code{\link[gear]{evgram}} function.
#' 
#' @param x An \code{evgram} object produced by the \code{\link[gear]{evgram}} function.
#' @param ... Not currently implemented.
#' @return NULL
#' @author Joshua French
#' @export
#' @examples 
#' data(co)
#' evgram(Al ~ 1, co, ~ easting + northing)
print.evgram = function(x, ...) {
  print(x$semivariogram)
}
