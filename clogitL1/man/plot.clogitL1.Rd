\name{plot.clogitL1}
\alias{plot.clogitL1}
\title{Plotting after fitting conditional logistic regression with elastic net penalties}
\description{
 Takes a \code{clogitL1} object and plots the parameter profile associated with it.
}


\usage{
 \method{plot}{clogitL1} (x, logX=T,
	add.legend=F, add.labels=T,
 	lty=1:ncol(x$beta), col=1:ncol(x$beta), ...) 
}

\arguments{
\item{x}{an object of type \code{clogitL1}.}
\item{logX}{should the horizontal axis be on log scale?}
\item{add.legend}{set to TRUE if legend should be printed in top right hand corner. Legend will contain names of variables in data.frame, if specified, otherwise will be numbered from 1 to p in order encountered in original input matrix x}
\item{add.labels}{set to TRUE if labels are to be added to curves at leftmost side. If variable names are available, these are plotted, otherwise, curves are numbered from 1 to p in order encountered in original input matrix x}
\item{lty}{usual 'lty' plotting parameter.}
\item{col}{usual 'col' plotting parameter.}
\item{...}{additional arguments to \code{plot} function}
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
plot(clObj, logX=TRUE)

# cross validation
clcvObj = cv.clogitL1(clObj)
plot(clcvObj)
}