\name{residuals.robmixglm}
\alias{residuals.robmixglm}
\title{Extract Model Residuals}
\usage{
\method{residuals}{robmixglm}(object,  type  =  c("deviance",  "pearson"),  \dots)
}
\arguments{
  \item{object}{an object for which the extraction of model residuals is
    meaningful.}
  \item{type}{Type of residual where valid types are deviance and pearson.}
  \item{\dots}{other arguments.}
}
\value{
  Residuals extracted from the object \code{object}.
}
\description{
  Extracts model residuals from objects returned by modeling functions.
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
BIC(forbes.robustmix)
plot(fitted(forbes.robustmix), residuals(forbes.robustmix))
}}
\keyword{models}
