\name{AIC}
\alias{AIC.robmixglm}
\title{AIC for robmixglm object}
\description{
Returns AIC for a robmixglm object.
}
\usage{
\method{AIC}{robmixglm}(object,  \ldots,  k  =  2)
}
\arguments{
  \item{object}{robmixglm object}
  \item{\dots}{ additional argument; currently none is used. }
  \item{k}{penalty per parameter}
}
\value{
AIC
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
AIC(forbes.robustmix)
}
}
\author{Ken Beath}
\keyword{methods}