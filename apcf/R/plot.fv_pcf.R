#' Plot PCF
#'
#' Plot method for the class "fv_pcf". Draws a pair correlation function and
#' a pointwise critical envelope if available.
#'
#' @param x an object of class [fv_pcf].
#' @param xlim,ylim the x and y limits of the plot. NULL indicates that the
#'        range of the finite values to be plotted should be used.
#' @param xticks,yticks points at which tick-marks are to be drawn.
#'        By default (when NULL) tickmark locations are computed
#' @param xlab,ylab a label for the x and y axis, respectively.
#' @param main,sub a main and sub title for the plot, see also [title()].
#' @param xaxis,yaxis a logical value or a 1-character string giving the
#'        desired type of axis. The following values are possible: "n" or FALSE
#'        for no axis, "t" for ticks only, "f" or TRUE for full axis, "o" for
#'        full axis in the outer margin.
#' @param ann a logical value indicating whether the default annotation
#'        (title and x and y axis labels) should appear on the plot.
#' @param bty a character string which determines the type of box which is
#'        drawn about the plotting region, see [par()].
#' @param ... additional parameter, currently without effect
#'
#' @return An object of class [fv_pcf] invisibly.
#'
#' @seealso [pat2dists()], [dists2pcf()]
#'
#' @examples
#' # it's advised against setting n_sim < 199
#' ds <- pat2dists(area=system.file("shapes/sim_area.shp", package="apcf"),
#'                 pattern=system.file("shapes/sim_pat_reg.shp", package="apcf"),
#'                 max_dist=25, n_sim=3)
#'
#' # derive PCF and envelope
#' pcf <- dists2pcf(ds, r=0.2, r_max=25, stoyan=0.15, n_rank=1)
#'
#' # a simple plot
#' plot(x=pcf, xlim=c(0, 20), ylim=c(0, 2.2))
#'
#' # a panel of four plots
#' op <- par(mfrow=c(2,2), oma=c(3,3,0,0), mar=c(0,0,2,2),
#'           mgp=c(2,0.5,0), tcl=-0.3)
#' plot(pcf, xaxis='t', yaxis='o', ann=FALSE)
#' plot(pcf, xaxis='t', yaxis='t', ann=FALSE)
#' plot(pcf, xaxis='o', yaxis='o', ann=FALSE)
#' plot(pcf, xaxis='o', yaxis='t')
#' par(op)
#'
#' @export
plot.fv_pcf <- function(x, xlim=NULL, ylim=NULL, xticks=NULL, yticks=NULL,
                        xlab="r", ylab="g(r)", main=NULL, sub=NULL,
                        xaxis=TRUE, yaxis=TRUE, ann=graphics::par("ann"),
                        bty='l', ...){
  stopifnot(is.fv_pcf(x))
  drawEnvelope <- "upr" %in% names(x)

  xlim <- if(is.null(xlim))
    range(0, x$r[is.finite(x$r)])
  else xlim
  ylim <- if(is.null(ylim))
    range(x$g[is.finite(x$g)])
  else ylim

  # set up plot
  graphics::plot.new()
  graphics::plot.window(xlim, ylim, xaxs="i", yaxs="i")

  xticks <- if(is.null(xticks))
    graphics::axTicks(1)
  else xticks
  yticks <- if(is.null(yticks))
    graphics::axTicks(2)
  else yticks

  # vertical lines
  graphics::segments(x0=xticks[-1], y0=ylim[1], x1=xticks[-1], y1=ylim[2],
                     col="gray", lty="11")

  # envelope with reference line at g(r)=1
  if(drawEnvelope){
    graphics::polygon(x=c(x$r, rev(x$r)), y=c(x$lwr, rev(x$upr)),
                      col="gray", border=NA)
    graphics::abline(h=1, col="white", lty="solid", lwd=1)
  } else {
    graphics::abline(h=1, col="black", lty="dotted", lwd=1)
  }

  # PCF
  graphics::lines(x$r, x$g, lwd=1.5)

  # axes and labels
  switch(ifelse(is.logical(xaxis),ifelse(xaxis, 'f', 'n'), xaxis),
         n = {},
         f = graphics::axis(1, at=xticks, lwd=0, lwd.ticks=1),
         t = graphics::axis(1, at=xticks, lwd=0, lwd.ticks=1, labels=FALSE),
         o = graphics::axis(1, at=xticks, lwd=0, lwd.ticks=1, outer=TRUE),
         stop("xaxis must be on of TRUE, FALSE, n, f, o, t")
  )
  switch(ifelse(is.logical(yaxis),ifelse(yaxis, 'f', 'n'), yaxis),
         n = {},
         f = graphics::axis(2, at=yticks, lwd=0, lwd.ticks=1),
         t = graphics::axis(2, at=yticks, lwd=0, lwd.ticks=1, labels=FALSE),
         o = graphics::axis(2, at=yticks, lwd=0, lwd.ticks=1, outer=TRUE),
         stop("yaxis must be on of TRUE, FALSE, n, f, o, t")
  )
  graphics::box(bty=bty)
  if(ann)
    graphics::title(main=main, sub=sub, xlab=xlab, ylab=ylab,
                    outer=if(xaxis == 'o' || yaxis == 'o'){TRUE} else {FALSE})

  invisible(x)
}
