#' Extract fitted values from a \code{geolm} object
#'
#' Extract the fitted values, i.e., the estimated mean
#' values, for an \code{object} produced by the
#' \code{\link[gear]{geolm}} function for a specified set of
#' covariates, \code{x}. If \code{x} is \code{NULL}, then
#' then \code{x} is taken from \code{object}.
#'
#' If the \code{object} has a known mean, i.e.,
#' \code{object$mu} is not \code{NULL}, then the function
#' returns the vector \code{rep(object$mu, m)}. If
#' \code{object} has estimated coefficients, then \code{x
#' \%*\% object$coeff} is returned.
#'
#' If \code{x} is missing, then \code{object$x} is used for
#' \code{x}. Naturally, \code{ncol(x)} must equal
#' \code{length(object$coeff)}. If \code{x} is \code{NULL}
#' and \code{object$mu} is not \code{NULL}, then \code{m} is
#' taken to be 1.
#'
#' @param object An object produced by the
#'   \code{\link[gear]{geolm}} function.
#' @param x A \eqn{m \times p} matrix of covariates for the
#'   locations where fitted values are desired. If
#'   \code{NULL}, \code{object$x} is used.
#' @param ... Not currently implemented.
#' @return The vector of fitted values.
#' @author Joshua French
#' @importFrom stats fitted
#' @export
#' @seealso \code{\link[stats]{fitted}}
#' @examples
#' data = data.frame(y = rnorm(10), x1 = runif(10),
#'                  x2 = runif(10))
#' d = as.matrix(dist(data[,c("x1", "x2")]))
#' mod = cmod_man(v = exp(-d), evar = 1)
#' gearmod = geolm(y ~ x1, data = data,
#'                 coordnames = ~ x1 + x2, mod = mod)
#' # fitted values for original observations
#' fitted(gearmod)
#' # fitted values for new observations
#' fitted(gearmod, x = cbind(1, rnorm(20)))
#' @rdname fitted.geolm
fitted.geolm = function(object, x, ...) {
  if (missing(x)) x = object$x
  arg_check_fitted_geolm(object, x)
  m = ifelse(is.null(x), 1, nrow(x))
  if (is.null(object$mu)) {
    return(c(x %*% object$coeff))
  } else {
    return(rep(object$mu, m))
  }
}

#' Check arguments of fitted.geolm
#'
#' @param object geolm object
#' @param x Matrix of covariates or NULL
#' @noRd
arg_check_fitted_geolm = function(object, x){
  if (is.null(object$mu)) {
    if (!is.matrix(x)) {
      stop("x must be a matrix")
    }
    if (nrow(x) == 0) {
      stop("x must have at least one row")
    }
    p = length(object$coeff)
    if (is.matrix(object$coeff)) { p = nrow(object$coeff) }
    if (ncol(x) != p) {
      stop("ncol(x) != number of regression coefficients")
    }
  }
}


