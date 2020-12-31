\name{fk_density}
\alias{fk_density}
\title{Fast univariate kernel density estimation}
\description{
  Uses recursive formulation for kernel sums as described in Hofmeyr (2019) to evaluate kernel estimate of the density exactly. Binning approximation also available for faster computation if needed. Default is exact evaluation on a grid, but evaluation at an arbitrary collection of points is possible.
}
\usage{
fk_density(x, h = 'Silverman', h_adjust = 1, beta = NULL, from = NULL,
              to = NULL, ngrid = 1000, nbin = NULL, x_eval = NULL)
}
\arguments{
  \item{x}{vector of sample points.}
  \item{h}{(optional) bandwidth to be used in estimate. Can be either positive numeric, or one of "Silverman" for Silverman's rule of thumb (Silverman, 1986) or "mlcv" for maximum pseudo-likelihood cross validation bandwidth. Default is Silverman's heuristic. Exact "mlcv" will be time consuming for samples of more than millions of points. Binning approximation for "mlcv" bandwidth needs at least 10000 bins for reasonable accuracy, and even more if density has very sharp features.}
  \item{h_adjust}{(optional) positive numeric. Final bandwidth will be h*h_adjust. Default value is 1. R's base function density uses Silverman's heuristic with h_adjust approximately 0.85.}
  \item{beta}{(optional) numeric vector of kernel coefficients. See Hofmeyr (2019) for details. The
  default is the smooth order one kernel described in the paper.}
  \item{from}{(optional) lower end of evaluation interval if evaluation on a grid is desired. Default is min(x)-6*h}
  \item{to}{(optional) upper end of evaluation interval if evaluation on a grid is desired. Default is max(x)+6*h}
  \item{ngrid}{(optional) integer number of grid points for evaluation. Default is 1000.}
  \item{nbin}{(optional) integer number of bins if binning estimator is to be used. The default is to compute the exact density on a grid of 1000 points.}
  \item{x_eval}{(optional) vector of evaluation points. The default if both ngrid and nbin are set to NULL is evaluation at the sample points themselves. If another specific set of points is required then ngrid must be set to null and x_eval supplied. Evaluation at arbitrary x_eval using binned approximation is also possible, in which case nbin and x_eval must both be supplied.}
}
\value{
  A named list with fields
  \item{$x}{the vector of points at which the density is estimated.}
  \item{$y}{the estimated density values.}
  \item{$h}{the value of the bandwidth.}
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.

Silverman, B. (1986) \emph{Density estimation for statistics and data analysis}, volume 26. CRC press.
}
\examples{
op <- par(no.readonly = TRUE)

set.seed(1)

### generate a bimodal Gaussian mixture sample with 100000 points

n1 <- rbinom(1, 100000, .5)

x <- c(rnorm(n1), rnorm(100000-n1)/4+2)



# --------- Example 1: Grid evaluation --------------#
### evaluate exact and binned approximation on a
### grid and plot along with true density. We can
### use both Silverman's heuristic and the "mlcv" estimate.

xs <- seq(-5, 3.5, length = 1000)

ftrue <- dnorm(xs)/2 + dnorm(xs, 2, 1/4)/2

fhat <- fk_density(x, from = -5, to = 3.5)

fhat_bin <- fk_density(x, nbin = 1000, from = -5, to = 3.5)

par(mfrow = c(2, 2))
plot(xs, ftrue, type = 'l', lwd = 4, col = rgb(.7, .7, .7),
    xlab = 'x', ylab = 'f(x)',
    main = 'Exact evaluation, Silverman bandwidth')
lines(fhat, lty = 2)

plot(xs, ftrue, type = 'l', lwd = 4, col = rgb(.7, .7, .7),
    xlab = 'x', ylab = 'f(x)',
    main = 'Binned approximation, Silverman bandwidth')
lines(fhat_bin, lty = 2)




fhat <- fk_density(x, from = -5, to = 3.5, h = 'mlcv')

fhat_bin <- fk_density(x, nbin = 1000, from = -5,
    to = 3.5, h = 'mlcv')

plot(xs, ftrue, type = 'l', lwd = 4, col = rgb(.7, .7, .7),
    xlab = 'x', ylab = 'f(x)',
    main = 'Exact evaluation, MLCV bandwidth')
lines(fhat, lty = 2)

plot(xs, ftrue, type = 'l', lwd = 4, col = rgb(.7, .7, .7),
    xlab = 'x', ylab = 'f(x)',
    main = 'Binned approximation, MLCV bandwidth')
lines(fhat_bin, lty = 2)


par(op)


# --------- Example 2: Evaluation at sample --------------#
### evaluate exact and binned approximation at the sample.
### Note that the output will be in the order of the original
### sample, and also the number of points will be large. It is
### not advisable, therefore, to simply plot these. We instead
### compute the mean squared deviations from the true density

ftrue <- sapply(x, function(xi) dnorm(xi)/2 + dnorm(xi, 2, 1/4)/2)

fhat <- fk_density(x, ngrid = NULL)

fhat_bin <- fk_density(x, nbin = 1000, x_eval = x)


mean((ftrue-fhat$y)^2)

mean((ftrue-fhat_bin$y)^2)


### now for MLCV bandwidth

fhat <- fk_density(x, h = 'mlcv', ngrid = NULL)

fhat_bin <- fk_density(x, nbin = 1000,
    h = 'mlcv', x_eval = x)


mean((ftrue-fhat$y)^2)

mean((ftrue-fhat_bin$y)^2)

}
\keyword{file}
