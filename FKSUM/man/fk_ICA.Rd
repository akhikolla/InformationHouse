\name{fk_ICA}
\alias{fk_ICA}
\title{Independent component analysis with sample entropy estimated via kernel density}
\description{
  Performs minimum entropy projection pursuit
}
\usage{
fk_ICA(X, ncomp = 1, beta = c(.25, .25), hmult = 1.5, it = 20, nbin = NULL)
}
\arguments{
  \item{X}{numeric data matrix (num_data x num_dimensions).}
  \item{ncomp}{integer number of indpenedent components to extract.}
  \item{beta}{numeric vector of kernel coefficients. See Hofmeyr (2019) for details. The
  default is the smooth order one kernel described in the paper.}
  \item{hmult}{positive numeric. The bandwidth in the kernel density is set to hmult multiplied by Silverman's rule of thumb value, which is based on the AMISE minimiser when the underlying distribution is Gaussian.}
  \item{it}{integer maximum number of iterations. The default is 20.}
  \item{nbin}{integer number of bins if binning estimator is to be used. The default is to compute the exact entropy estimate from the kde.}
}
\value{
  A named list with fields
  \item{$X}{the data matrix given as argument.}
  \item{$K}{the pre-whitening matrix.}
  \item{$W}{the estimated un-mixing matrix.}
  \item{$S}{the estimated source matrix, S = sweep(X,2,colMeans(X),'-')\%*\%K\%*\%W.}
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
op <- par(no.readonly = TRUE)

set.seed(1)

### Generate a set of data with 5 non-Gaussian components
### and 20 noise (Gaussian) components. The non-Gaussian
### components are cubed Gaussian, exponential, uniform,
### t3 and bimodal Gaussian mixture.

errdim = 10
ndat = 2000
X = cbind(rnorm(ndat)^3/5, rexp(ndat), runif(ndat)*sqrt(12)*2, rt(ndat, 3),
    c(rnorm(ndat/2), rexp(ndat/2)+2), matrix(rnorm(errdim*ndat), ncol = errdim))

### Generate a random mixing matrix and mixed data XX

R = matrix(rnorm((5+errdim)*(5+errdim)), ncol = 5+errdim)
XX = X\%*\%R

### Apply fk_ICA to XX to extract components.
### Note that ordering of extracted components is
### in some sense arbitrary w.r.t. their generation

model <- fk_ICA(XX, ncol(XX))

par(mfrow = c(1, 5))

for(i in 1:5) plot(density(model$S[,i]))

par(op)
}
\keyword{file}
