#' Plot \code{geardf} object
#'
#' Plot a \code{geardf} object produced by the
#' \code{\link[gear]{geardf}} or \code{predict.geolm*}
#' functions. See the \code{\link[autoimage]{autopoints}}
#' and \code{\link[autoimage]{autoimage}} functions for
#' advanced options.
#'
#' @param x A \code{geardf} object produced by the
#'   \code{\link[gear]{geardf}} or \code{predict.geolm*}
#'   functions.
#' @param zcol The names of the columns of \code{x} to plot.
#'   The coordinate columns are automatically stripped from
#'   \code{zcol}.
#' @param interp A logical value indicating whether the
#'   values should be interpolated onto a grid. If
#'   \code{FALSE}, then the
#'   \code{\link[autoimage]{autopoints}} function is used to
#'   construct a heated scatterplot. If \code{TRUE}, then
#'   the \code{\link[autoimage]{autoimage}} function is used
#'   to create a heat map. The default is \code{FALSE}.
#' @param common.legend A logical value indicating whether a
#'   common legend should be used for the
#'   scatterplots/images. The default is \code{FALSE}.
#' @param ... Additional arguments to passed to the
#'   \code{\link[autoimage]{autopoints}} or
#'   \code{\link[autoimage]{autoimage}} functions.
#' @return NULL
#' @author Joshua French
#' @export
#' @method plot geardf
#' @seealso \code{\link[autoimage]{autopoints}},
#'   \code{\link[autoimage]{autoimage}}
#' @examples
#' data(toydata)
#' # newdata must have columns with prediction coordinates
#' newdata = data.frame(x1 = runif(10), x2 = runif(10))
#'
#' # specify a standard covariance model
#' mod = cmod_std(model = "exponential", psill = 1, r = 1)
#'
#' # geolm for universal kriging
#' geolm_uk = geolm(y ~ x1 + x2, data = toydata, mod = mod,
#'              coordnames = c("x1", "x2"))
#'
#' # prediction for universal kriging
#' pred_uk = predict(geolm_uk, newdata, return_type = "geardf")
#' # heated scatterplot
#' plot(pred_uk)
#' # interpolated image of results
#' plot(pred_uk, interp = TRUE)
#' # plot only predictions and rmspe with different colors
#' plot(pred_uk, c("pred", "rmspe"), col = cm.colors(5))
#' #'plot only predictions with coarser interpolation grid
#' plot(pred_uk, "pred", interp = TRUE,
#'      interp.args = list(no.X = 10, no.Y = 10))
plot.geardf = function(x, zcol = names(x), interp = FALSE,
                       common.legend = FALSE, ...) {
  if (length(interp) != 1 | !is.logical(interp)) {
    stop("interp must be a single logical value")
  }
  coordnames = attr(x, "coordnames")
  znames = setdiff(zcol, coordnames)
  arglist = list(...)
  arglist$common.legend = common.legend
  if (is.null(arglist$xlab)) {
    arglist$xlab = coordnames[1]
  }
  if (is.null(arglist$ylab)) {
    arglist$ylab = coordnames[2]
  }
  if (is.null(arglist$main)) {
    arglist$main = znames
  }
  arglist$x = x[, coordnames[1], drop = TRUE]
  arglist$y = x[, coordnames[2], drop = TRUE]
  arglist$z = x[, znames, drop = FALSE]

  if (interp) {
    do.call(autoimage::autoimage, arglist)
  } else {
    do.call(autoimage::autopoints, arglist)
  }
}
