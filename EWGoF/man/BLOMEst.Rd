\name{BLOMEst}
\alias{BLOMEst}
\title{
Blom's estimators of the two parameters of the Weibull distribution 
}
\description{
Computes the Blom's estimators of the shape and scale parameters of the Weibull distribution from an i.i.d sample x. 
It also gives the sample \eqn{\check{y}} after using the logarithmic transformation (\eqn{\check{y}=(\check{shape})\ln(x/\check{scale})}, where \eqn{\check{shape}} and \eqn{\check{scale}} are the estimated shape and scale parameters). 
}
\usage{
BLOMEst(x)
}
\arguments{
  \item{x}{a numeric vector of data values.}
}
\value{
A list containing the following elements:
  \item{eta}{the Blom's estimator of the scale parameter of the Weibull distribution (\eqn{\check{scale}}).}
  \item{beta}{the Blom's estimator of the shape parameter of the Weibull distribution (\eqn{\check{shape}}).}
  \item{y}{the pseudo-observations \eqn{\check{y}} after using the logarithmic transformation and the Blom's estimators. }
}
\references{
Blom G., Statistical Estimates and Transformed Beta-variables. New York: Wiley, 1958.}
\author{Meryam KRIT}

\examples{
x <- rweibull(50,2,3)

#Value of the Blom's estimator of the scale parameter
BLOMEst(x)$eta

#Value of the Blom's estimator of the shape parameter
BLOMEst(x)$beta

}
\keyword{Blom's estimators}
\keyword{linear unbiased estimators}
