\name{h_K_to_Gauss}
\alias{h_K_to_Gauss}
\title{Bandwidth conversion to Gaussian}
\description{
  A naive and simple, but useful way of transforming a bandwidth for use with a kernel implemented in FKSUM to an appropriate value to be
  used with the Gaussian kernel. The transformation if based on the relative values for the kernels for the AMISE minimal
  bandwidth for the purpose of density estimation. The transformation has been useful for other problems as well.
}
\usage{
h_K_to_Gauss(h, beta)
}
\arguments{
  \item{h}{positive numeric bandwidth value appropriate for use with the FKSUM kernel with coefficients beta to an appropriate value
  for use with the Gaussian kernel.}
  \item{beta}{numeric vector of kernel coefficients.}
}
\value{
  positive numeric value, the bandwidth to be used for estimation using the Gaussian kernel.
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
### Use the package for data driven bandwidth for use with
### the Gaussian kernel in existing implementations

### generate sample from bimodal Gaussian mixture with varying scale

x <- c(rnorm(100000), rnorm(100000)/4 + 2)

### estimate bandwidth with the package using MLCV and convert

bwml <- h_K_to_Gauss(fk_density(x, h = 'mlcv',
    beta = c(.25,.25))$h, c(.25,.25))

bwml_binned <- h_K_to_Gauss(fk_density(x, h = 'mlcv',
    beta = c(.25,.25), nbin = 10000)$h, c(.25,.25))

xs <- seq(-3, 4, length = 1000)
plot(xs, dnorm(xs)/2+dnorm(xs,2,1/4)/2, type = 'l',
    lwd = 4, col = rgb(.7, .7, .7))

lines(density(x, bw = bwml), lty = 2, lwd = 2)

lines(density(x, bw = bwml_binned), lty = 3, lwd = 2)
}
\keyword{file}
