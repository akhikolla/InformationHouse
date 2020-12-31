\name{MEst}
\alias{MEst}
\title{
Moment Estimators of the two parameters of the Weibull distribution 
}
\description{
This function computes the Moment Estimators (MEs) of the shape and scale parameters of the Weibull distribution from an i.i.d sample x. It also gives the sample \eqn{\breve{y}} after using the logarithmic transformation (\eqn{\breve{y}=(\breve{shape})\ln(x/\breve{scale})}, where \eqn{\breve{shape}} and \eqn{\breve{scale}} are the estimated shape and scale parameters).}
\usage{
MEst(x)
}
\arguments{
  \item{x}{a numeric vector of data values.}
}
\details{
The elements of the numeric vector should be positive. The support of the Weibull distribution is R+*. These estimators are based on the moments of the extreme value distribution.
}
\value{
A list containing the following elements:
  \item{eta}{the moment estimator of the scale parameter of the Weibull distribution (\eqn{\breve{scale}}).}
  \item{beta}{the moment estimator of the shape parameter of the Weibull distribution (\eqn{\breve{shape}}).}
  \item{y}{the pseudo-observations \eqn{\breve{y}} after using the logarithmic transformation and the MEs.}
}
\references{
Rinne H., \emph{ The Weibull distribution - A handbook}, CRC-Chapman & Hall, 2009. 

}
\author{Meryam KRIT}
\examples{
x <- rweibull(50,2,3)

Est <- MEst(x)

#Value of the moment estimator of the scale parameter
Est$eta

#Value of the moment estimator of the shape parameter
Est$beta
}

\keyword{Moment estimators}
