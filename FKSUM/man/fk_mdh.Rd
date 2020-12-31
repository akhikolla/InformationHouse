\name{fk_mdh}
\alias{fk_mdh}
\title{Minimum density hyperplanes}
\description{
  Estimates minimum density hyperplanes for clustering using projection pursuit, based on the method of Pavlidis et al. (2016)
}
\usage{
fk_mdh(X, v0 = NULL, hmult = 1, beta = c(.25,.25), alphamax = 1)
}
\arguments{
  \item{X}{numeric data matrix (num_data x num_dimensions).}
  \item{v0}{(optional) vector representing the initial projection direction. Default is the first principal direction.}
  \item{hmult}{(optional) positive numeric. The bandwidth in the kernel density is set to hmult multiplied by Silverman's rule of thumb value,
  which is based on the AMISE minimiser when the underlying distribution is Gaussian.}
  \item{beta}{(optional) numeric vector of kernel coefficients. The default is the smooth order one kernel described by Hofmeyr (2019).}
  \item{alphamax}{(optional) maximum/final (scaled) distance of the optimal hyperplane from the mean of the data.
  The default is alphamax = 1.}
}
\value{
  A named list with fields
  \item{$v}{the optimal projection vector.}
  \item{$b}{the location of the minimum density hyperplane orthogonal to v.}
}
\references{
Pavlidis N.G., Hofmeyr D.P., Tasoulis S.K. (2016) "Minimum Density Hyperplanes", \emph{Journal of Machine Learning Research}, 17(156), 1–33.

Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
op <- par(no.readonly = TRUE)

set.seed(1)

### Generate data from a simple 10 component mixture model in 10 dimensions:
### Determine means of components

mu <- matrix(runif(100), ncol = 10)

### Determine scales of components (diagonal elements of covariance)

covs <- matrix(rexp(100), ncol = 10)/10

### Determine cluster indicator matrix

I <- t(rmultinom(2000, 1, 1:10))

### Determine mean and residual matrix

M <- I\%*\%mu

R <- matrix(rnorm(20000), 2000, 10)*(I\%*\%covs)

### Data is given by the sum of these

X <- M + R

### Find the minimum density hyperplane separator and plot
### the projected data as well as the PCA plot for comparison

mdh <- fk_mdh(X)

par(mfrow = c(2, 2))

plot(X\%*\%mdh$v, X\%*\%eigen(cov(X))$vectors[,2], xlab = 'mdh vector',
    ylab = '2nd principal component', main = 'MDH solution')
abline(v = mdh$b, col = 2)
plot(X\%*\%eigen(cov(X))$vectors[,1:2], xlab = '1st principal component',
    ylab = '2nd principal component', main = 'PCA plot')

plot(fk_density(X\%*\%mdh$v))
abline(v = mdh$b, col = 2)
plot(fk_density(X\%*\%eigen(cov(X))$vectors[,1]))

par(op)
}
\keyword{file}
