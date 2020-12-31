\name{plot.unpencoxIC}
\alias{plot.unpencoxIC}
\title{Plot method for \code{unpencoxIC} object}
\usage{
\method{plot}{unpencoxIC}(x, what = "cum.hazard", xlim, ylim, xlab, ylab, axes = FALSE, ...)
}
\arguments{
\item{...}{for S4 method only.}
\item{x}{An object of class \code{unpencoxIC} returned by the \code{unpencoxIC} function.}
\item{what}{A character string specifying which function will be plotted. Default is "cum.hazard", which plots the estimated baseline cumulative hazard function. Set to "survival" to plot the estimated baseline survival function.}
\item{xlim}{A vector with two elements for the limits of follow-up time.}
\item{ylim}{A vector with two elements for the limits of y-axis.}
\item{xlab}{A label for the x axis.}
\item{ylab}{A label for the y axis.}
\item{axes}{A logical value drawing both axes. Default is \code{FALSE}.}
}
\description{
Plot method for \code{unpencoxIC} object for plotting the estimated baseline culmulative function and the estimated baseline survival function.
}
\details{
The \code{x} argument must be the object returned by the \code{unpencoxIC} function. Note that \code{plot} provides the conditional survival function for left truncated data, which is analogous to the function (5) of Alioum and Commenges (1996). See the usages in the examples given below.
}

\examples{
  library(ALassoSurvIC)

  data(ex_ICLT) # interval censored and left truncated data
  lowerIC <- ex_ICLT$lowerIC
  upperIC <- ex_ICLT$upperIC
  trunc <- ex_ICLT$trunc
  X <- ex_ICLT[, -c(1:3)]
  result <- unpencoxIC(lowerIC, upperIC, X, trunc)

  plot(result)  # plotting the estimated baseline cumulative hazard function by default
  plot(result, what = "survival")  # plotting the estimated baseline survival function

}

\references{Alioum, A. and Commenges, D. (1996). A proportional hazards model for arbitrarily censored and truncated data. Biometrics 52, 512-524.
}

\keyword{methods}
