#' Print a \code{geolm} object
#'
#' Print an \code{object} produced by
#'  \code{\link[gear]{geolm}}.
#'
#' @param x An object produced by the
#'   \code{\link[gear]{geolm}} function.
#' @param ... Not currently implemented.
#' @return Information about the \code{geolm}
#' @author Joshua French
#' @export
#' @seealso \code{\link[base]{print}}
#' @noRd
#' @examples
#' data = data.frame(y = rnorm(10), x1 = runif(10),
#'                  x2 = runif(10))
#' d = as.matrix(dist(data[,c("x1", "x2")]))
#' mod = cmod_man(v = exp(-d), evar = 1)
#' gearmod = geolm(y ~ x1, data = data, mod = mod,
#'                 coordnames = ~ x1 + x2)
#' gearmod
print.geolm = function(x, ...) {
  # cat("Call:\n")
  # print(x$call)
  cat("Call: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (is.null(x$mu)) {
    cat("Coefficients:\n")
    print(x$coeff)
    cat("\n")
  } else {
    cat(paste("Mu: ", x$mu, "\n\n"))
  }
  cat(paste("Mod class:", class(x$mod)[1]))
}
