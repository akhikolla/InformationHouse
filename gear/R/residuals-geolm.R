#' Extract residuals from a \code{geolm} object
#'
#' Extract the residuals for an \code{object} produced by
#'  the \code{\link[gear]{geolm}}.
#'
#' @param object An object produced by the
#'   \code{\link[gear]{geolm}} function.
#' @param ... Not currently implemented.
#' @return The vector of residuals.
#' @author Joshua French
#' @importFrom stats residuals
#' @export
#' @seealso \code{\link[stats]{residuals}}
#' @examples
#' data = data.frame(y = rnorm(10), x1 = runif(10),
#'                  x2 = runif(10))
#' d = as.matrix(dist(data[,c("x1", "x2")]))
#' mod = cmod_man(v = exp(-d), evar = 1)
#' gearmod = geolm(y ~ x1, data = data,
#'                 coordnames = ~ x1 + x2, mod = mod)
#' # fitted values for original observations
#' residuals(gearmod)
#' @rdname residuals.geolm
residuals.geolm = function(object, ...) {
  c(object$resid)
}
