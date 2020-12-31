#' Update linear model for geostatistical data
#'
#' \code{update} updates a geostatistical linear model based
#' on the given model.
#'
#' @param object An object produced by the \code{geolm}
#'   function.
#' @param mod A spatial dependence model object obtained
#'   from one of the \code{cmod_*} functions or the
#'   \code{estimate} function.
#' @param ... Not implemented.
#' @return Returns an object of the same class as
#'   \code{object}.
#'
#' @author Joshua French
#' @seealso \code{\link[stats]{update}}
#' @examples
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#'
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' coords = cbind(x1, x2)
#' psill = 2 # partial sill
#' r = 4 # range parameter
#' evar = .1 # error variance
#' fvar = .1 # add finescale variance
#' # one can't generally distinguish between evar and fvar, but
#' # this is done for illustration purposes
#'
#' cmod_std = cmod_std("exponential", psill = psill, r = r,
#'                     evar = evar, fvar = fvar)
#'
#' cmod_std2 = cmod_std("exponential", psill = psill + 1,
#'                      r = r + .5, evar = evar + .01,
#'                      fvar = fvar)
#'
#' # check geolm update for universal kriging
#' gear1 = geolm(y ~ x1 + x2, data = data, mod = cmod_std,
#'               coordnames = c("x1", "x2"))
#'
#' gear2 = geolm(y ~ x1 + x2, data = data, mod = cmod_std2,
#'               coordnames = c("x1", "x2"))
#' gear2b = update(gear1, cmod_std2)
#' gear2$call = NULL
#' gear2b$call = NULL
#' identical(gear2, gear2b)
#' @rdname update.geolm
#' @export
update.geolm = function(object, mod, ...) {
  message("Update.geolm is internally governed by update.geolm_cmodMan, update.geolm_cmodStd, etc.")
}
