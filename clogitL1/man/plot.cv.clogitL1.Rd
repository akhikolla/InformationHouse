\name{plot.cv.clogitL1}
\alias{plot.cv.clogitL1}
\title{Plotting after cross validating conditional logistic regression with elastic net penalties}
\description{
 Takes a \code{cv.clogitL1} object and plots the CV deviance curve with standard error bands and minima.
}

\usage{
 \method{plot}{cv.clogitL1} (x, ...)
}

\arguments{
\item{x}{an object of type \code{cv.clogitL1}.}
\item{...}{additional arguments to \code{plot} function}
}


\references{
\url{http://www.jstatsoft.org/v58/i12/}
}

\seealso{
 \code{\link{cv.clogitL1}}
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
plot(clcvObj)
}