\name{MLEst}
\alias{MLEst}
\title{
Maximum Likelihood Estimators of the two parameters of the Weibull distribution 
}
\description{
This function computes the Maximum Likelihood Estimators (MLEs) of the shape and scale parameters of the Weibull distribution from an i.i.d sample x. 
It also gives the sample \eqn{\hat{y}} after using the logarithmic transformation (\eqn{\hat{y}=(\hat{shape})\ln(x/\hat{scale})}, where \eqn{\hat{shape}} and \eqn{\hat{scale}} are the estimated shape and scale parameters).
}
\usage{
MLEst(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values.}
}
\details{
The elements of the numeric vector should be positive. The support of the Weibull distribution is R+*.
}
\value{
A list containing the following elements:
  \item{eta}{the maximum likelihood estimator of the scale parameter of the Weibull distribution (\eqn{\hat{scale}}).}
  \item{beta}{the maximum likelihood estimator of the shape parameter of the Weibull distribution (\eqn{\hat{shape}}).}
  \item{y}{the pseudo-observations \eqn{\hat{y}} after using the logarithmic transformation and the MLEs. }
}
\references{
D'Agostino R.B. and Stephens M.A., \emph{Goodness-of-fit techniques}, Marcel Dekker, 1986.

Krit M., Gaudoin O., Xie M. and Remy E., Simplified likelihood goodness-of-fit tests for the Weibull distribution, \emph{Communications in Statistics - Simulation and Computation.}
}
\author{Meryam KRIT}

\examples{
x <- rweibull(50,2,3)

Est <- MLEst(x)

#Value of the maximum likelihood estimator of the scale parameter
Est$eta

#Value of the maximum likelihood estimator of the shape parameter
Est$beta
}

\keyword{Maximum likelihood estimators}
