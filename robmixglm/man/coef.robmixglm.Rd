\name{coef}
\alias{coef}
\alias{coef.robmixglm}
\title{Coefficients for a robmixglm object}
\description{
Returns coefficients for a robmixglm object. Only the coefficients for the linear part of the model are returned. Additional coefficients may be obtained using summary().
}
\usage{
\method{coef}{robmixglm}(object,  \ldots)
}
\arguments{
  \item{object}{robmixglm object}
  \item{\dots}{additional argument; currently none is used. }
}
\value{
coef
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
coef(forbes.robustmix)
}}
\author{Ken Beath}
\keyword{methods}