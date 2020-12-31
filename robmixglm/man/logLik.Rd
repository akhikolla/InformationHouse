\name{logLik}
\alias{logLik.robmixglm}
\title{log Likelikelihood for robmixglm object}
\description{
Returns log Likelihood for a robmixglm object.
}
\usage{
\method{logLik}{robmixglm}(object,  \ldots)
}
\arguments{
  \item{object}{robmixglm object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
The loglikelihood.
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
logLik(forbes.robustmix)
}}
\author{Ken Beath}
\keyword{methods}