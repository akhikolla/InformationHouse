#' Contrast adjustments for RasterLayer plots
#'
#' A plotting function that provides methods for improving the contrast
#' between values.
#'
#'
#' @param x a RasterLayer
#' @param type character. Possible values are "linear", "hist.equal", and "sd"
#' @param n number of standard deviations to include if \code{type = "sd"}
#' @param ... Additional arguments for raster::plot
#'
#' @details
#' If \code{type = "hist.equal"}, a histogram equalization procedure will be
#' applied to the values of \code{x}. If \code{type = "sd"}, the values of
#' \code{x} will be scaled between values that fall between \code{n} standard
#' deviations of the mean.
#'
#' @examples
#' mod <- enfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' sm <- sensitivity_map(mod)
#' stretchPlot(sm)
#' stretchPlot(sm, type = "hist.equal")
#' stretchPlot(sm, type = "sd", n = 2)
#'
#' @importFrom stats ecdf sd
# @importFrom raster calc cellStats
# @importMethodsFrom raster [ plot
#'
#' @export

setGeneric("stretchPlot", function(x, type = "linear", n, ...) {
  standardGeneric("stretchPlot")
})

#' @rdname stretchPlot
setMethod("stretchPlot",
          signature(x = "RasterLayer"),
          function(x, type = "linear", n, ...) {

            xx <- pretty(values(x), 2)
            p <- length(xx)

            if(type == "hist.equal") {
              y <- .stretch(x, type = "hist.equal")
              plot(y[[1]], breaks = .quantile_breaks(y[[1]]),
                   axis.args = list(at = y[[2]][-p],
                                    labels = xx[-p]), ...)
            }
            if(type == "linear") {
              y <- .stretch(x, type = "linear")
              plot(y[[1]],
                   axis.args = list(at = y[[2]],
                                    labels = xx), ...)
            }
            if(type == "sd") {
              y <- .stretch(x, type = "sd", n = n)
              f <- duplicated(y[[2]])
              plot(y[[1]],
                   axis.args = list(at = y[[2]][!f],
                                    labels = xx[!f]), ...)
            }
          }
)


#' @keywords internal
.stretch <- function(x, type = "hist.equal", n) {
  if (type == "hist.equal") {
    ecdfun <- stats::ecdf(getValues(x))
    y <- calc(x, fun = function(x) round(ecdfun(x)*255, 0))
    v <- pretty(values(x), 2)
    vv <- round(ecdfun(v)*255, 0)
    return(list(y, vv))
  } else if (type == "sd") {
    if (missing(n)) stop("number of standard deviations not specified")
    x.sd <- cellStats(x, sd)
    x.mean <- cellStats(x, mean)
    x.max <- x.mean + x.sd*n
    x.min <- x.mean - x.sd*n
    sdfun <- function(y) {
      y[y < 0] <- x.min
      y[y > x.max] <- x.max
      tt <- (y - x.min) / (x.max - x.min)
      return(tt*255)
    }
    y <- calc(x, fun = function(x) sdfun(x))
    v <- pretty(values(x), 2)
    vv <- sdfun(v)
    return(list(y, vv))
  } else if (type == "linear") {
    y <- x
    vv <- pretty(values(x), 2)
    return(list(y, vv))
  }

}

#' @keywords internal
.quantile_breaks <- function(ras){
  qts <- quantile(ras)
  seq1 <- seq(qts[1], qts[2], length.out = 65)
  seq2 <- seq(qts[2], qts[3], length.out = 65)
  seq3 <- seq(qts[3], qts[4], length.out = 65)
  seq4 <- seq(qts[4], qts[5], length.out = 64)
  brks <- c(seq1, seq2[-1], seq3[-1], seq4[-1])
  brks
}
