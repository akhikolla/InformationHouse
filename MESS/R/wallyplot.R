#' @rdname wallyplot
#' @export
wallyplot.default <- function(x, y=x, FUN=residualplot,
                              hide=TRUE,
                              simulateFunction=rnorm,
                              ...) {

  simulateFunction <- match.fun(simulateFunction)

  if (is.vector(y) && length(x)==length(y))
    y <- cbind(y, sapply(1:8, function(k) {simulateFunction(length(x))}))

  if (!is.numeric(x) || !is.matrix(y))
    stop("x and y input must be a vector and matrix, respectively")

  if (length(x) != nrow(y))
    stop("x and y does not conform")

  if (ncol(y) != 9)
    stop("y must be a matrix with 9 columns")

  if (!is.numeric(x) && !is.numeric(y)) {
    stop("Must have a pair of numeric vectors or an lm object as input")
  }

  cc <- complete.cases(x, y)
  x <- x[cc]
  y <- y[cc,]

  plot3x3(x,y, FUN=FUN, hide=hide, ...)

  invisible(NULL)
}


#' @import graphics
#' @rdname wallyplot
#' @export
wallyplot.lm <- function(x, y=x, FUN=residualplot,
	                 hide=TRUE,
                         simulateFunction=rnorm,
                         ...) {

  # Extract information from model fit
  y <- rstandard(x)
  x <- predict(x)

  wallyplot.default(x, y, FUN=FUN, hide=hide, simulateFunction=simulateFunction, ...)
}




#' Plots a Wally plot
#'
#' Produces a 3x3 grid of residual- or qq-plots plots from a lm object. One of
#' the nine subfigures is the true residual plot/qqplot while the remaining are
#' plots that fulfill the assumptions of the linear model
#'
#' Users who look at residual plots or qqnorm plots for the first time often
#' feel they lack the experience to determine if the residual plot is okay or
#' if the model assumptions are indeed violated. One way to convey "experience"
#' is to plot a series of graphical model validation plots simulated under the
#' model assumption together with the corresponding plot from the real data and
#' see if the user can pinpoint one of them that looks like an odd-one-out. If
#' the proper plot from the real data does not stand out then the assumptions
#' are not likely to be violated.
#'
#' The Wallyplot produces a 3x3 grid of plots from a lm object or from a set of
#' pairs of x and y values. One of the nine subfigures is the true plot while
#' the remaining are plots that fulfill the assumptions of the linear model.
#' After the user interactively hits a key the correct residual plot
#' (correponding to the provided data) is shown.
#'
#' The plotting function can be set using the \code{FUN} argument which should
#' be a function that accepts \code{x}, \code{y} and \code{...} arguments and
#' plots the desired figure. When \code{y} is a single vector the same length
#' as \code{x} then the function \code{simulateFunction} is used to generate
#' the remaining y values corresponding the situations under the null.
#'
#' For a description of the features of the default residual plot see the help page for \code{\link{residualplot}}.
#'
#' @aliases wallyplot wallyplot.lm wallyplot.default
#' @param x a numeric vector of x values, or an lm object.
#' @param y a numeric vector of y values of the same length as x or a n * 9
#' matrix of y values - one column for each of the nine plots to make. The
#' first column is the one corresponding to the results from the dataset
#' @param FUN a function that accepts an \code{x}, \code{y} and \code{...}
#' argument and produces a graphical model validation plots from the \code{x}
#' and \code{y} values.
#' @param hide logical; if \code{TRUE} (the default) then the identity of the
#' true residual plot is hidden until the user presses a key. If \code{FALSE}
#' then the true residual plot is shown in the center.
#' @param simulateFunction The function used to produce y values under the null
#' hypothesis. Defaults to rnorm
#' @param ... Other arguments passed to the plot function \code{FUN}
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2014) \emph{Teaching 'Instant Experience' with
#' Graphical Model Validation Techniques}. Teaching Statistics (36), p 23-26
#' @keywords iplot
#' @examples
#'
#' \dontrun{
#' data(trees)
#' res <- lm(Volume ~ Height + Girth, data=trees)
#' wallyplot(res)
#'
#'
#' # Create a grid of QQ-plot figures
#' # Define function to plot a qq plot with an identity line
#' qqnorm.wally <- function(x, y, ...) { qqnorm(y, ...) ; abline(a=0, b=1) }
#' wallyplot(res, FUN=qqnorm.wally, main="")
#'
#' # Define function to simulate components+residuals for Girth
#' cprsimulate <- function(n) {rnorm(n)+trees$Girth}
#' # Create the cpr plotting function
#' cprplot <- function(x, y, ...) {plot(x, y, pch=20, ...) ;
#'                                  lines(lowess(x, y), lty=3)}
#' # Create the Wallyplot
#' wallyplot(trees$Girth, trees$Girth+rstudent(res), FUN=cprplot,
#'           simulateFunction=cprsimulate, xlab="Girth")
#' }
#'
#' @export
wallyplot <- function(x, y=x, FUN=residualplot,
                      hide=TRUE,
                      simulateFunction=rnorm,
                       ...) {
  UseMethod("wallyplot")
}


# qqnorm.wally <- function(x, y, ...) {qqnorm(y, ...) ; abline(a=0, b=1)}


plot3x3 <- function(x, y, FUN=plot, hide=TRUE, ylim=range(y), mar=c(4, 4, .1, .1)+.1, ...) {
  # Input check

  if (!is.numeric(x) || !is.matrix(y))
    stop("x and y input must be a vector and matrix, respectively")

  if (length(x) != nrow(y))
    stop("x and y does not conform")

  if (ncol(y) != 9)
    stop("y must be a matrix with 9 columns")

  oldpar <- par(no.readonly = TRUE)
  par(mar=mar)

  FUN <- match.fun(FUN)

  pos <- c(1:9)
  if (hide)
    pos <- sample(pos)

  par(mfrow=c(3,3))

  for (i in 1:length(pos)) {
    FUN(x, y[,pos[i]], ylim=ylim, ...)
  }

  if (hide) {
    readline("Hit <Enter> to show the original plot. ")
  }

  figpos <- order(pos)[1]
  par(mfg=c(figpos %/% 3.1 + 1, figpos - (figpos %/% 3.1)*3 ))
  box(col="red", lwd=2)

  par(oldpar)
}
