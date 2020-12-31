\name{outlierTest}
\alias{outlierTest}
\alias{outlierTest.robmixglm}
\title{
Test for the presence of outliers.}
\description{
Uses the parametric bootstrap to test for the presence of outliers.}
\usage{
outlierTest(object,  R = 999,  cores  =  max(detectCores() - 1,  1))
}
\arguments{
  \item{object}{
A robmixglm object with a mixture (robust) random effects distribution.
}
  \item{R}{number of bootstrap replications}
  \item{cores}{Number of cores to be used in parallel. Default is one less than available.}
}
\value{
An outlierTest object which is the object of class \dQuote{boot} returned by the call to \command{boot}.
}
\details{
Performs a parametric bootstrap to compare models with and without outliers.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
hospcosts.robustmix <- robmixglm(costs~adm+age+dest+ins+loglos+sex, family = "gamma", 
    data = hospcosts, cores = 1)
summary(hospcosts.robustmix)
summary(outlierTest(hospcosts.robustmix,  cores  =  1))
}}
\keyword{method}
