\name{summary.cv.clogitL1}
\alias{summary.cv.clogitL1}
\title{Summary after cross validation of conditional logistic regression with elastic net penalties}
\description{
 Provides summary of conditional logistic regression models after cross validation
}

\usage{
 \method{summary}{cv.clogitL1} (object, ...)
}

\arguments{
\item{object}{an object of type \code{cv.clogitL1} for which the summary is to be produced.}
\item{...}{additional arguments to \code{summary} method.}
}

\details{
Extracts pertinent information from the supplied \code{cv.clogitL1} objects. See below for details on output value.
}

\value{
A list with the following fields:

\item{lambda_minCV}{value of regularisation parameter minimising CV deviance}
\item{beta_minCV}{coefficient profile at the minimising value of the regularisation parameter. Whole dataset used to compute estimates.}
\item{nz_beta_minCV}{number of non-zero coefficients in the CV deviance minimising coefficient profile.}
\item{lambda_minCV1se}{value of regularisaion parameter minimising CV deviance (using 1 standard error rule)}
\item{beta_minCV1se}{coefficient profile at the 1-standard-error-rule value of the regularisation parameter. Whole dataset used to compute estimates.}
\item{nz_beta_minCV1se}{number of non-zero coefficients in the 1-standard-error-rule coefficient profile.}

}

\references{
\url{http://www.jstatsoft.org/v58/i12/}
}

\seealso{
 \code{\link{clogitL1}}, \code{\link{plot.cv.clogitL1}}
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
plot(clObj, logX=TRUE)

# cross validation
clcvObj = cv.clogitL1(clObj)
summary(clcvObj)
}