\name{fitted.robmixglm}
\alias{fitted.robmixglm}
\title{
Fitted values.}
\description{
Calculates the fitted values.}
\usage{
\method{fitted}{robmixglm}(object, \ldots)
}
\arguments{
  \item{object}{
A robmixglm object with a mixture (robust) random effects distribution.}
  \item{\ldots}{Other parameters. (not used)}
}
\value{
A vector of the fitted values.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
BIC(forbes.robustmix)
plot(fitted(forbes.robustmix), residuals(forbes.robustmix))
}}
\keyword{method}
