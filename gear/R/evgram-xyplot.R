#' Plot \code{evgram} object
#'
#' Plots \code{evgram} object produced by the
#' \code{\link[gear]{evgram}} function using the
#' \code{\link[lattice]{xyplot}} function. Note: the
#' \code{lattice} package must be loaded (i.e.,
#' \code{library(lattice)} or \code{lattice::xyplot}
#' must be specifically called for this function to work.
#' See Examples.
#'
#' @param x An \code{evgram} object produced by the \code{\link[gear]{evgram}} function.
#' @param ... Additional arguments to pass the \code{\link[lattice]{xyplot}} function to change aspects of the plot.
#' @param split A logical value indicating whether, for a directional
#' semivariogram, the directional semivariograms should be displayed
#' in a single or split panels.  Default is FALSE, for a single panel.
#' @return NULL
#' @author Joshua French
#' @export
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link[gear]{evgram}}
#' @examples
#' data(co)
#' v = evgram(Al ~ 1, co, ~ easting + northing)
#' if (requireNamespace("lattice")) {
#'    lattice::xyplot(v)
#' }
#' v2 = evgram(Al ~ 1, co, ~ easting + northing, angle = 22.5, ndir = 4)
#' if (requireNamespace("lattice")) {
#'    lattice::xyplot(v2)
#'    # show how attributes can be changed using different
#'    # arguments available in lattice::xyplot.
#'    lattice::xyplot(v2, col = 2:5)
#'    lattice::xyplot(v2, col = 2:5, pch = 1:4)
#'    lattice::xyplot(v2, col = 2:5, pch = 1:4, lty = 2:5, type = "b")
#'    lattice::xyplot(v2, col = 2:5, pch = 1:4, lty = 2:5, type = "b",
#'     key=list(text=list(levels(as.factor(v2$semi$angle))),
#'     space='right', points=list(pch=1:4, col=2:5),
#'     lines=list(col=2:5, lty = 2:5)))
#'    lattice::xyplot(v2, split = TRUE)
#' }
xyplot.evgram = function(x, ..., split = FALSE) {
  if (!requireNamespace("lattice")) {
    stop("lattice must be installed to enable this functionality")
  }
  arg_check_split(split)
  largs = list(...)
  if (is.null(largs$xlim)) {
    largs$xlim = c(0, max(x$semivariogram$distance))
  }
  if (is.null(largs$lim)) {
    largs$ylim = c(0, max(x$semivariogram$semivariance))
  }
  largs$data = x$semivariogram
  if (x$ndir == 1)   {
    largs$x = stats::formula("semivariance ~ distance")
  }
  else if (split) {
    largs$x = stats::formula("semivariance ~ distance | angle")
  } else {
    largs$x = stats::formula("semivariance ~ distance")
    largs$groups = x$semivariogram$angle
  }
  do.call(lattice::xyplot, largs)
}
