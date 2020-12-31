\name{outlierProbs}
\alias{outlierProbs}
\alias{outlierProbs.robmixglm}
\title{
Calculate outlier probabilities for each observation.}
\description{
For the normal mixture random effect calculates the probability that each observation is an outlier based on the posterior probability of it being an outlier.}
\usage{
outlierProbs(object)
}
\arguments{
  \item{object}{
A metaplus object with a mixture (robust) random effects distribution.
}
}
\value{
\item{outlier.prob}{Posterior probability that each observation is an outlier}
}
\details{
The outlier probabilities are obtained as the posterior probabilities of each observation being an outlier based on the fitted mixture model.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
outlierProbs(forbes.robustmix)
}}
\keyword{method}
