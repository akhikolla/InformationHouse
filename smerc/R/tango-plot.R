#' Plots an object of class \code{tango}.
#'
#' Plots results of \code{\link{tango.test}}.  If Monte
#' Carlo simulation was not used to produce \code{x}, then a  a density plot of
#' the (approximate) null distribution of \code{tstat.chisq} is produced, along
#' with a vertical line for the observed \code{tstat}.
#' If a Monte Carlo test was used to produce \code{x}, then a scatterplot of
#' the \code{gof.sim} versus \code{sa.sim} is compared to the observed values
#' \code{gof} and \code{sa}, respectively.
#'
#' @param x An object of class \code{tango} to be plotted.
#' @param ... Additional graphical parameters passed to \code{plot} function.
#' @param obs.list A list containing arguments for the
#' \code{\link[graphics]{points}} function, which is used to
#' plot the \code{gof} and \code{sa} components,
#' when appropriate.
#' @param sim.list A list containing arguments for the
#' \code{\link[graphics]{points}} function, which is used to
#' plot the \code{gof.sim} and \code{sa.sim} components,
#' when appropriate.
#' @method plot tango
#' @export
#' @seealso \code{\link{tango.test}}
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("x", "y")])
#' w = dweights(coords, kappa = 1)
#' x1 = tango.test(nydf$cases, nydf$pop, w)
#' plot(x1)
#' x2 = tango.test(nydf$cases, nydf$pop, w, nsim = 49)
#' plot(x2)
plot.tango = function(x, ...,
                      obs.list = list(pch = 20),
                      sim.list = list(pch = 2)) {
  if (class(x) != "tango") {
    stop("x should be a tango object")
  }

  if (!is.null(x$gof.sim)) {
    gof = c(x$gof.sim, x$gof)
    sa = c(x$sa.sim, x$sa)
    graphics::plot(sa, gof, type = "n", ...)
    sim.list$x = x$sa.sim
    sim.list$y = x$gof.sim
    do.call(graphics::points, sim.list)
    obs.list$x = x$sa
    obs.list$y = x$gof
    do.call(graphics::points, obs.list)
  } else {
    tstat = seq(0, max(x$tstat.chisq, stats::qchisq(.99, df = x$dfc)))
    density = stats::dchisq(tstat, df = x$dfc)
    graphics::plot(tstat, density, type = "l", ...)
    graphics::abline(v = x$tstat.chisq)
  }
}
