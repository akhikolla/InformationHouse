\name{print.outlierTest}
\alias{print.outlierTest}
\alias{summary.outlierTest}
\alias{print.summary.outlierTest}
\title{Print an outlierTest object}
\description{
Print an outlierTest object.
}
\usage{
\method{print}{outlierTest}(x,  \ldots)
}
\arguments{
  \item{x}{outlierTest object}
  \item{\ldots}{further arguments (not currently used)}
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
summary(forbes.robustmix)
print(outlierTest(forbes.robustmix,  cores  =  1))
}
}
\author{Ken Beath}
\keyword{methods}