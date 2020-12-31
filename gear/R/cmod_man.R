#' Manual covariance models for geostatistical data.
#'
#' \code{cmod_man} manually creates a covariance model
#' object (\code{cmodMan}) for geostatistical data.
#'
#' Note that \code{v} includes the error variance, i.e.,
#' \code{v = vz + ve}, where \code{vz} is the covariance
#' matrix of the filtered process and the variance matrix
#' of the errors is \code{ve = diag(evar/weights)}, where
#' the \code{weights} come from the \code{geolm} object the
#' \code{cmodMan} object is associated with.
#'
#' @param v The covariance matrix of the observed data,
#'   including any errors.  The matrix should be square,
#'   symmetric, and positive definite, though that latter
#'   two conditions are not checked.
#' @param evar The variance of the errors.  Must be
#'   non-negative number.  The default is 0.
#'
#' @return Returns a \code{cmodMan} object.
#'
#' @author Joshua French
#' @export
#' @examples
#' coords = matrix(runif(6), ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod_man(v = exp(-d), evar = 1)
cmod_man = function(v, evar = 0) {
  arg_check_cmod_man(v, evar)
  return(structure(list(v = v, evar = evar), class = "cmodMan"))
}

#' Argument check cmod_man
#'
#' @param v The covariance matrix of the observed data, including any errors.  The matrix should be square, symmetric, and positive definite, though that latter two conditions are not checked.
#' @param evar The variance of the errors.  Must be non-negative number.  The default is 0.
#' @return NULL
#' @noRd
arg_check_cmod_man = function(v, evar) {
  if (!is.numeric(v) | nrow(v) != ncol(v) ) stop("v should be a square matrix")
  arg_check_evar(evar)
}


