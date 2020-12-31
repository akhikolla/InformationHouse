#' Hanging rootogram for normal distribution
#' 
#' Create a hanging rootogram for a quantitative numeric vector and compare it
#' to a Gaussian distribution.
#' 
#' The mean and standard deviation of the Gaussian distribution are calculated
#' from the observed data unless the \code{mu} and \code{s} arguments are
#' given. %% ~~ If necessary, more details than the description above ~~
#' 
#' @aliases rootonorm rootogram
#' @param x a numeric vector of values for which the rootogram is desired
#' @param breaks Either the character string \sQuote{Sturges} to use Sturges'
#' algorithm to decide the number of breaks or a positive integer that sets the
#' number of breaks.
#' @param type if \code{"hanging"} then a hanging rootogram is plotted, and if
#' \code{"deviation"} then deviations from zero are plotted.
#' @param scale The type of transformation. Defaults to \code{"sqrt"} which
#' takes square roots of the frequencies. \code{"raw"} yields untransformed
#' frequencies.
#' @param zeroline logical; if \code{TRUE} a horizontal line is added at zero.
#' @param linecol The color of the density line for the normal distribution.
#' The default is to make a \code{red} density line.
#' @param rectcol a colour to be used to fill the bars.  The default of
#' \code{lightgray} yields lightgray bars.
#' @param xlab,ylab plot labels.  The \code{xlab} and \code{ylab} refer to the
#' x and y axes respectively
#' @param yaxt Should y axis text be printed. Defaults to \code{n}.
#' @param ylim the range of y values with sensible defaults.
#' @param mu the mean of the Gaussian distribution. Defaults to the sample mean
#' of \code{x}.
#' @param s the standard deivation of the Gaussian distribution. Defaults to
#' the sample std.dev. of \code{x}.
#' @param gap The distance between the rectangles in the histogram.
#' @param \dots further arguments and graphical parameters passed to
#' \code{plot}.
#' @return Returns a vector of counts of each bar. This may be changed in the
#' future. The plot is the primary output of the function.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Tukey, J. W. 1972. \emph{Some Graphic and Semigraphic Displays}.
#' In \emph{Statistical Papers in Honor of George W. Snedecor}, p. 293-316.
#' @keywords hplot
#' @examples
#' 
#' oldpar <- par()
#' par(mfrow=c(2,2))
#' rootonorm(rnorm(200))
#' rootonorm(rnorm(200), type="deviation", scale="raw")
#' rootonorm(rnorm(200), mu=1)
#' rootonorm(rexp(200), mu=1)
#' par(oldpar)
#' 
#' @export rootonorm
rootonorm <- function(x, breaks="Sturges",
                      type=c("hanging", "deviation"),
                      scale=c("sqrt", "raw"),
                      zeroline=TRUE,                      
                      linecol="red", rectcol="lightgrey",
                      xlab=xname,
                      ylab="Sqrt(frequency)",
                      yaxt="n",                      
                      ylim=NULL,
                      mu=mean(x), s=sd(x),
                      gap=0.1, ...) {


  if (!is.numeric(x)) 
    stop("'x' must be numeric")

  # Fix the xlabel if it isn't specified
  xname <- deparse(substitute(x))
 
  scale <- match.arg(scale)
  if (is.character(scale) && scale == "raw") {
    scale <- match.fun("as.numeric")
    if (missing(ylab)) {
      ylab <- "Frequency"
    }
  } else {
    scale <- match.fun(scale)
  } 

  type <- match.arg(type)
  
  h <- hist(x, breaks=breaks, plot=FALSE)
  if (!h$equidist) stop("breaks must be equally spaced")
  
  nbins <- length(h$counts)
  nobs <- sum(h$counts)

  expected <- nobs*diff(pnorm(h$breaks, mu, s))
  
  d.gap <- min(diff(h$breaks)) * gap /2
 
  plot.range <- range(pretty(h$breaks))  
  z <- seq(plot.range[1], plot.range[2], length.out=200)
  z.y <- min(diff(h$breaks))*nobs*dnorm(z, mu, s)
  
  minval <- min(scale(expected)-scale(h$counts))

  if (is.null(ylim)) {
    ylim <- c(minval, scale(max(expected,z.y)))
  }
  
  plot(z, z, type="n",
       xlab=xlab,
       ylab=ylab,
       yaxt=yaxt,
       ylim=ylim,
       ...)

  if (type=="deviation") {  
    for(i in 1:nbins) {
      rect(h$breaks[i]+d.gap, scale(expected[i])-scale(h$counts[i]),
           h$breaks[i+1]-d.gap, 0, col=rectcol)
    }    
  } else {
    for(i in 1:nbins) {
      rect(h$breaks[i]+d.gap, scale(expected[i])-scale(h$counts[i]),
           h$breaks[i+1]-d.gap, scale(expected[i]), col=rectcol)
    }
  }
  
  lines(z, scale(z.y), col=linecol, ...)
  if (zeroline) {
    abline(h=0, lty=3)
  }

  invisible(h$counts)
}
