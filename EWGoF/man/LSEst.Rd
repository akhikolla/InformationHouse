\name{LSEst}
\alias{LSEst}
\title{
Least Squares Estimators of the two parameters of the Weibull distribution
}
\description{
This function computes the Least Squares Estimators (LSEs) of the shape and scale parameters of the Weibull distribution, based on the probability plot, from an i.i.d sample x. 
It also gives the sample \eqn{\tilde{y}} after using the logarithmic transformation (\eqn{\tilde{y}=(\widetilde{shape})\ln(x/\widetilde{scale})}, where \eqn{\widetilde{shape}} and \eqn{\widetilde{scale}} are the estimated shape and scale parameters). 
}
\usage{
LSEst(x)
}
\arguments{
  \item{x}{a numeric vector of data values.}
}
\details{
The elements of the numeric vector should be positive. The support of the Weibull distribution is R+*. These estimators are used by Liao and Shimokawa; they are based on the probability plot and symmetrical ranks.
}
\value{
A list containing the following elements:
  \item{eta}{the least squares estimator of the scale parameter of the Weibull distribution (\eqn{\widetilde{scale}}).}
  \item{beta}{the least squares estimator of the shape parameter of the Weibull distribution (\eqn{\widetilde{shape}}).}
  \item{y}{the pseudo-observations \eqn{\tilde{y}} after using the logarithmic transformation and the LSEs.}
}
\references{
Liao M. and Shimokawa T., A new goodness-of-fit test for type-I extreme-value and 2-parameter Weibull distributions with estimated parameters, \emph{Journal of Statistical Computation and Simulation}, 64 (1), 23-48, 1999.

Krit M., Gaudoin O., Xie M. and Remy E., Simplified likelihood goodness-of-fit tests for the Weibull distribution, \emph{Communications in Statistics - Simulation and Computation}.
}
\author{Meryam KRIT}

\examples{
x <- rweibull(50,2,3)

Est <- LSEst(x)

#Value of the least squares estimator of the scale parameter
Est$eta

#Value of the lest squares estimator of the shape parameter
Est$beta
}
\keyword{Liao-shimokawa}
\keyword{Least squares estimators}
