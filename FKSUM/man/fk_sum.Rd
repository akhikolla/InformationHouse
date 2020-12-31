\name{fk_sum}
\alias{fk_sum}
\title{Fast Exact Kernel Sum Evaluation}
\description{
  Computes exact (and binned approximations of) kernel and kernel derivative sums with arbitrary weights/coefficients. Computation is based on the method of Hofmeyr (2019).
}
\usage{
fk_sum(x, omega, h, x_eval = NULL, beta = c(.25,.25),
    nbin = NULL, type = "ksum")
}
\arguments{
  \item{x}{numeric vector of sample points.}
  \item{omega}{numeric vector of weights.}
  \item{h}{numeric bandwidth (must be strictly positive).}
  \item{x_eval}{vector of evaluation points. Default is to evaluate at the sample points themselves.}
  \item{beta}{numeric vector of kernel coefficients. Default is c(.25, .25); the smooth order 1 kernel.}
  \item{nbin}{integer number of bins for binning approximation. Default is to compute the exact sums.}
  \item{type}{one of "ksum": returns the kernel sums, "dksum": returns the kernel derivative sums and "both": returns a matrix cbind(ksum, dksum).}
}
\value{
  A vector (if type\%in\%c("ksum","dksum")) of kernel sums, or kernel derivative sums. A matrix (if type == "both") with kernel sums in its first column and kernel derivative sums in its second column.
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
### Compute density estimates directly with
### kernel sums and constant normalising weights

set.seed(1)
n <- 150000
num_Gauss <- rbinom(1, n, 2 / 3)
x <- c(rnorm(num_Gauss), rexp(n - num_Gauss) + 1)

hs <- seq(.025, .1, length = 5)
xeval <- seq(-4, 8, length = 1000)
ftrue <- 2 / 3 * dnorm(xeval) + 1 / 3 * dexp(xeval - 1)
plot(xeval, ftrue, lwd = 6, col = rgb(.8, .8, .8), xlab = "x",
  ylab = "f(x)", type = "l")

for(i in 1:5) lines(xeval, fk_sum(x, rep(1 / hs[i] / n, n), hs[i],
    x_eval = xeval), lty = i)
}
\keyword{file}
