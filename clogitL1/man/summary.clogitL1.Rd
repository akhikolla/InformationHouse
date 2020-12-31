\name{summary.clogitL1}
\alias{summary.clogitL1}
\title{Summary after fitting conditional logistic regression with elastic net penalties}
\description{
 Takes a \code{clogitL1} object and produces a summary of the sequence of models fitted.
}

\usage{
 \method{summary}{clogitL1} (object, ...)
}

\arguments{
\item{object}{an object of type \code{clogitL1}.}
\item{...}{any additional arguments passed to \code{summary} method}
}

\details{
Returns a list with a elements \code{Coefficients}, which holds the matrix of coefficients estimated (each row holding the estimates for a given value of the smoothing parameter) and \code{Lambda}, which holds the vector of smoothing parameters at which fits were produced.
}

\references{
\url{http://www.jstatsoft.org/v58/i12/}
}

\seealso{
 \code{\link{clogitL1}}
}

\examples{
set.seed(145)

# data parameters
K = 10 # number of strata
n = 5 # number in strata
m = 2 # cases per stratum
p = 20 # predictors

# generate data
y = rep(c(rep(1, m), rep(0, n-m)), K)
X = matrix (rnorm(K*n*p, 0, 1), ncol = p) # pure noise
strata = sort(rep(1:K, n))

par(mfrow = c(1,2))
# fit the conditional logistic model
clObj = clogitL1(y=y, x=X, strata)
summary(clObj)
}