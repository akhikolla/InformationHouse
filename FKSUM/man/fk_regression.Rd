\name{fk_regression}
\alias{fk_regression}
\title{Fast univariate kernel regression}
\description{
  Uses recursive formulation for kernel sums as described in Hofmeyr (2019) for exact evaluation of kernel regression estimate.  Both local linear and Nadaraya-Watson (Watson, 1964; Nadaraya, 1964) estimators are available. Binning approximation also available for faster computation if needed. Default is exact evaluation on a grid, but evaluation at sample points is also possible.
}
\usage{
fk_regression(x, y, h = 'amise', beta = NULL, from = NULL,
    to = NULL, ngrid = 1000, nbin = NULL, type = 'loc-lin')
}
\arguments{
  \item{x}{vector of covariates.}
  \item{y}{vector of response values.}
  \item{h}{(optional) bandwidth to be used in estimate. Can be either positive numeric, or one of "amise" for a rough estimate of the asymptotic mean integrated square error minimiser, or "cv" for leave-one-out cross validation error minimiser based on squared error. Default is "amise". Cross validation not available for binning approximation.}
  \item{beta}{(optional) numeric vector of kernel coefficients. See Hofmeyr (2019) for details. The
  default is the smooth order one kernel described in the paper.}
  \item{from}{(optional) lower end of evaluation interval if evaluation on a grid is desired. Default is min(x)}
  \item{to}{(optional) upper end of evaluation interval if evaluation on a grid is desired. Default is max(x)}
  \item{ngrid}{(optional) integer number of grid points for evaluation. Default is 1000.}
  \item{nbin}{(optional) integer number of bins if binning estimator is to be used. The default is to compute the exact density on a grid of 1000 points.}
  \item{type}{(optional) one of "loc-lin" and "NW" if local linear or Nadaraya-Watson is desired, respectively. Default is local linear estimator.}
}
\value{
  A named list with fields
  \item{$x}{the vector of points at which the regression function is estimated.}
  \item{$y}{the estimated function values.}
  \item{$h}{the value of the bandwidth.}
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.

Nadaraya, E.A. (1964) "On estimating regression." \emph{Theory of Probability \& Its Applications},9(1), 141–142.

Watson, G.S. (1964) "Smooth regression analysis." \emph{Sankhya: The Indian Journal of Statistics}, Series A, pp. 359–372.
}
\examples{
op <- par(no.readonly = TRUE)

set.seed(1)

n <- 2000
x <- rbeta(n, 2, 2) * 10
y <- 3 * sin(2 * x) + 10 * (x > 5) * (x - 5)
y <- y + rnorm(n) + (rgamma(n, 2, 2) - 1) * (abs(x - 5) + 3)

xs <- seq(min(x), max(x), length = 1000)
ftrue <- 3 * sin(2 * xs) + 10 * (xs > 5) * (xs - 5)

fhat_loc_lin <- fk_regression(x, y)
fhat_NW <- fk_regression(x, y, type = 'NW')

par(mfrow = c(2, 2))
plot(fhat_loc_lin, main = 'Local linear estimator with amise bandwidth')

plot(fhat_NW, main = 'NW estimator with amise bandwidth')


fhat_loc_lin <- fk_regression(x, y, h = 'cv')
fhat_NW <- fk_regression(x, y, type = 'NW', h = 'cv')

plot(fhat_loc_lin, main = 'Local linear estimator with cv bandwidth')

plot(fhat_NW, main = 'NW estimator with cv bandwidth')

par(op)
}
\keyword{file}
