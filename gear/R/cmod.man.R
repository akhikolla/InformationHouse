#' Standard covariance models for geostatistical data.
#'
#' \code{cmod.man} manually creates a covariance matrix
#' object (\code{cmodMan}) for geostatistical data. This
#' function will be deprecated in the future. Please update
#' your code to use the \code{\link[gear]{cmod_man}}
#' function.
#'
#' Note that \code{v} includes the error variance, i.e.,
#' \code{v = vz + ve}, where \code{vz} is the covariance
#' matrix of the filtered process and the variance matrix
#' of the errors is \code{ve = diag(evar/weights)}, where
#' the \code{weights} come from the \code{geolm} object the
#' \code{cmodMan} object is associated with.
#'
#' @inheritParams cmod_man
#'
#' @return Returns a \code{cmodMan} object.
#'
#' @author Joshua French
#' @export
#' @examples
#' coords = matrix(runif(20), ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod.man(v = exp(-d), evar = 1)
cmod.man = function(v, evar = 0) {
  warning("This function will be deprecated in the future. Please update your code to use the cmod_man function.")
  cmod_man(v = v, evar = evar)
}

