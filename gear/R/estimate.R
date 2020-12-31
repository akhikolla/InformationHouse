#' Estimate model parameters
#'
#' \code{estimate} estimates the parameters of the specified
#' model. The function is written to automatically adapt
#' based on the class of \code{object}.  Currently, the
#' \code{estimate} function performs maximum likelihood
#' estimation for objects produced by the
#' \code{\link[gear]{geolm}} function.
#'
#' @param object A model object produced by the
#'   \code{\link[gear]{geolm}} function.
#' @param ... Currently unimplemented
#' @author Joshua French
#' @export
#' @examples
#' set.seed(10)
#' n = 100
#' @seealso \code{\link[gear]{estimate.geolm_cmodStd}},
#' @export
estimate = function(object, ...) {
  UseMethod("estimate")
}
