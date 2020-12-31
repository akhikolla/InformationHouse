plot.unpencoxIC <- function(x, what = "cum.hazard", xlim, ylim, xlab, ylab, axes = FALSE, ...) {

  match.call()

  smallest.trunc <- x$smallest.trunc
  dotnames <- names(list(...))
  if (any(dotnames == "type")) stop("The graphical argument 'type' is not allowed.")
  if (missing(what)) {
    what <- "cum.hazard"
  } else {
    if(!(what %in% c("cum.hazard", "survival"))) stop("The 'what' argument must be 'cum.hazard' or 'survival'.")
  }

  lambda <- x$lambda
  nlambda <- length(lambda)
  lambda.set <- x$lambda.set

  x.coordinate <- c(smallest.trunc, c(lambda.set[, 2]))

  if (missing(xlab)) xlab <- "Time"


  if(what == "cum.hazard") {
    y.coordinate <- c(0, cumsum(lambda))
    if (missing(ylab)) ylab <- "Baseline Cumulative Hazard Function"
  }
  if(what == "survival") {
    y.coordinate <- c(1, exp(-cumsum(lambda)))
    if (missing(ylab)) ylab <- "Baseline Survival Function"
  }

  if(missing(xlim)) xlim <- c(smallest.trunc, max(x.coordinate))
  if(missing(ylim)) ylim <- c(0, max(y.coordinate))

  x.max <- xlim[2]
  y.max <- ylim[2]

  plot(x.coordinate, y.coordinate, type = "s", xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, axes = axes, ...)

  if(axes == FALSE) {
    x.axis0 <- signif(seq(from = smallest.trunc, to = x.max, length.out = 6),2)
    x.axis <- c(x.axis0[x.axis0 <= signif(x.max, 2)], signif(x.max, 2))
    y.axis0 <- signif(seq(from = 0, to = y.max, length.out = 6),2)
    y.axis <- c(y.axis0[y.axis0 <= signif(y.max, 2)], signif(y.max, 2))
    axis(1, x.axis)
    axis(2, y.axis)
  }

}
