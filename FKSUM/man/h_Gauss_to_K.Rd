\name{h_Gauss_to_K}
\alias{h_Gauss_to_K}
\title{Bandwidth conversion from Gaussian}
\description{
  A naive and simple, but useful way of transforming a bandwidth for use with the Gaussian kernel to an appropriate value to be used
  with a kernel implemented in the FKSUM package. The transformation if based on the relative values for the kernels for the AMISE minimal
  bandwidth for the purpose of density estimation. The transformation has been useful for other problems as well.
}
\usage{
h_Gauss_to_K(h, beta)
}
\arguments{
  \item{h}{positive numeric bandwidth value appropriate for use with the Gaussian kernel.}
  \item{beta}{numeric vector of kernel coefficients.}
}
\value{
  positive numeric value, the bandwidth to be used for estimation using the kernel with coefficients beta.
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
### Use existing bandwidth selection with the package

### generate sample from bimodal Gaussian mixture with varying scale

x <- c(rnorm(100000), rnorm(100000)/4 + 2)

### estimate bandwidth using bw.SJ

bwsj <- h_Gauss_to_K(bw.SJ(x), c(.25,.25))

xs <- seq(-3, 4, length = 1000)
plot(xs, dnorm(xs)/2+dnorm(xs,2,1/4)/2, type = 'l',
    col = rgb(.7, .7, .7), lwd = 4)

lines(fk_density(x, h = bwsj), lty = 2, lwd = 2)

}
\keyword{file}
