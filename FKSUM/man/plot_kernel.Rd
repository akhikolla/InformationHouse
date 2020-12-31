\name{plot_kernel}
\alias{plot_kernel}
\title{Plot the shape of a kernel function implemented in FKSUM based on its vector of beta coefficients}
\description{
  Plots the kernel function for a given set of beta coefficients. NB: coefficients will be normalised so that the kernel describes
  a probability density.
}
\usage{
plot_kernel(beta, type = 'l', ...)
}
\arguments{
  \item{beta}{positive numeric vector of kernel coefficients.}
  \item{type}{(optional) plot type, as in base plot() function. The default is a line plot.}
  \item{...}{(optional) any additional arguments accepted by base plot() function.}
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
### Plot order 4 smooth kernel

plot_kernel(1/factorial(0:4))

### Use a different line style

plot_kernel(1/factorial(0:4), lty = 2)

}
\keyword{file}
