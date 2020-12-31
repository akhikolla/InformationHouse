#' Fit a \code{geolm}
#'
#' \code{geolm_fit} fits a \code{geolm} based on the
#' specified \code{mod}. This is effectively an internal
#' function.
#'
#' @param x The matrix of covariates.
#' @param y The vector of observed responses.
#' @param coords The coordinates of the observed data set.
#' @param weights A vector that is proportional to the reciprocal variances of the errors.
#' @param coordnames A vector of length 2 with the names of the columns in \code{data} containing the coordinates, e.g., \code{c("long", "lat")}.
#' @param n The number of observations.
#' @param call The \code{\link[base]{match.call}} for the \code{geolm}.
#' @param coeff_names A character string with the variable names associated with the coefficents.
#' @inheritParams geolm
#'
#' @return Returns a \code{geolm} with necessary
#'   components needed for \code{\link[gear]{estimate}} and
#'   \code{\link{predict}}.
#' @author Joshua French
#' @export
#' @examples
#' # no examples since you shouldn't be using this function!
geolm_fit = function(mod, x, y, coords, mu, weights, formula, coordnames, n, call, coeff_names) {
  UseMethod("geolm_fit")
}
