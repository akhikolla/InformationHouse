#' Plot \code{evgram} object
#'
#' Plot an \code{evgram} object produced by the
#' \code{\link[gear]{evgram}} function. The plotting
#' function internally calls the
#' \code{\link[ggplot2]{autoplot}} function. Note: the
#' \code{ggplot2} package must be loaded (i.e.,
#' \code{library(autplot)} or \code{ggplot2::autoplot}
#' must be specifically called for this function to work.
#' See Examples.
#'
#' @param object An \code{evgram} object produced by the
#'   \code{\link[gear]{evgram}} function.
#' @param ... Not used
#' @param split A logical value indicating whether, for a
#'   directional semivariogram, the directional
#'   semivariograms should be displayed in a single or split
#'   panels.  Default is \code{FALSE}, for a single panel.
#' @return NULL
#' @author Joshua French
#' @export
#' @examples
#' data(co)
#' v = evgram(Al ~ 1, co, ~ easting + northing)
#' if (requireNamespace("ggplot2")) {
#'    ggplot2::autoplot(v)
#' }
#' v2 = evgram(Al ~ 1, co, ~ easting + northing, angle = 22.5, ndir = 4)
#' # ggplot2 must manually be loaded for this to work
#' if (requireNamespace("ggplot2")) {
#'    ggplot2::autoplot(v2)
#'    ggplot2::autoplot(v2, split = TRUE)
#' }
autoplot.evgram = function(object, ..., split = FALSE) {
  if (!requireNamespace("ggplot2")) {
    stop("ggplot2 must be installed to enable this functionality")
  }
  arg_check_split(split)
  # appease R CMD check
  distance = semivariance = angle = NULL
  # map aesthetics based on split
  # construct ggplot object
  g = ggplot2::ggplot(data = object$semivariogram) +
    ggplot2::xlim(0, max(object$semivariogram$distance)) +
    ggplot2::ylim(0, max(object$semivariogram$semivariance))
  if (object$ndir == 1) {
    g +
      ggplot2::geom_point(ggplot2::aes(x = distance,
                                       y = semivariance))
  } else if (split) {
    g +
      ggplot2::geom_point(ggplot2::aes(x = distance,
                                       y = semivariance)) +
      ggplot2::facet_wrap(~ angle)
  } else {
    g +
      ggplot2::geom_point(ggplot2::aes(x = distance,
                                       y = semivariance,
                                       pch = angle,
                                       col = angle))
  }
}
