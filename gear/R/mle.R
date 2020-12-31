#' Finds maximum likelihood estimates of model parameters
#' for a geostatistical model
#'
#' \code{mle} estimates the parameters of a geostatistical
#' linear model.
#'
#' @param object A geostatistical linear model object
#'   producted by the \code{geolm} function.
#' @param ... Currently unimplemented
#' @author Joshua French
#' @export
#' @examples
#' set.seed(10)
#' n = 100
#' @seealso \code{\link[gear]{estimate}},
#' @export
mle = function(object, ...) {
  warning("mle will be deprecated in a future version of this package. Please use the estimate function instead.")
  warning("You may need to change the 'est' argument to 'noise_type', as the name of this argument has changed. See the estimate function.")
  estimate(object, ...)
}
