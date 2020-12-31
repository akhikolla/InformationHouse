\name{BIC}
\alias{BIC.robmixglm}
\title{BIC for robmixglm object}
\description{
Returns BIC for a robmixglm object.
}
\usage{
\method{BIC}{robmixglm}(object,  \ldots)
}
\arguments{
  \item{object}{robmixglm object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
BIC
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
BIC(forbes.robustmix)
}}
\author{Ken Beath}
\keyword{methods}