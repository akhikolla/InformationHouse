\name{print.clogitL1}
\alias{print.clogitL1}
\title{Printing after fitting conditional logistic regression with elastic net penalties}
\description{
 Takes a \code{clogitL1} object and prints a summary of the sequence of models fitted.
}

\usage{
 \method{print}{clogitL1} (x, digits = 6, ...)
}

\arguments{
\item{x}{an object of type \code{clogitL1}.}
\item{digits}{the number of significant digits after the decimal to be printed}
\item{...}{additional arguments to \code{print} function}
}

\details{
prints a 3 column data frame with columns:
\itemize{
\item{\code{Df}:}{ number of non-zero parameters in model}
\item{\code{DevPerc}:}{ percentage of null deviance explained by current model}
\item{\code{Lambda}:}{ associated \eqn{\lambda} value}
}
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
clObj
}