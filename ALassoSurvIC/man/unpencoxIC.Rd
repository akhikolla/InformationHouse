\name{unpencoxIC}
\alias{unpencoxIC}
\alias{unpencoxIC.default}
\alias{print.unpencoxIC}
\title{Performing unpenalized nonparametric maximum likelihood estimation for interval censored and possibly left truncated data}
\usage{
\method{unpencoxIC}{default}(lowerIC, upperIC, X, trunc = NULL, normalize.X = TRUE,
  cl = NULL, tol = 0.001, niter = 1e+05, string.cen = Inf, string.missing = NA, ...)
}

\arguments{
\item{...}{for S4 method only.}
\item{lowerIC}{A numeric vector for the lower limit of the censoring interval.}

\item{upperIC}{A numeric vector for the upper limit of the censoring interval.}

\item{X}{A numeric matrix for the covariates that will be used for variable selection.}

\item{trunc}{A numeric vector for left truncated time. If supplied, the function performs the variable selection for interval censored and left truncated data. If \code{trunc} is missing, the data will be considered as interval censored data.}

\item{normalize.X}{A logical value: if \code{normalize.X} = \code{TRUE}, the covariate matrix \code{X} will be normalized before fitting models. Default is \code{TRUE}.}

\item{cl}{A cluster object created by \code{\link[parallel]{makeCluster}} in the \code{parallel} package. If \code{NULL}, no parallel computing is used by default. See details below.}

\item{tol}{A numeric value for the absolute iteration convergence tolerance.}

\item{niter}{A numeric value for the maximum number of iterations.}

\item{string.cen}{A string indicating right censoring for \code{upperIC}. Default is \code{Inf}.}

\item{string.missing}{A string indicating missing value. Default is \code{NA}.}


}
\description{
The \code{unpencoxIC} function performs unpenalized nonparametric maximum likelihood estimation. The function provides unpenalized nonparametric maximum likelihood estimates, standard errors, and 95\% confidence intervals. The full details are available in Li et al. (2019).
}

\details{
The cluster object, created by \code{\link[parallel]{makeCluster}} in the \code{parallel} package, can be supplied with the \code{cl} argument to reduce computation time via parallel computing. The parallel computing will be used when calculating the hessian matrix of the estimates. How to use the parallel computing is illustrated in one of the examples given below.

Use the \code{baseline} function and the \code{plot} function to extract and plot the estimate of the baseline cumulative hazard function, respectively, from the object returned by the \code{unpencoxIC}. The \code{plot} function also provides the plot for estimated baseline survival function. See the usages in the examples given below.
}

\references{Li, C., Pak, D., & Todem, D. (2019). Adaptive lasso for the Cox regression with interval censored and possibly left truncated data. Statistical methods in medical research. \href{https://doi.org/10.1177/0962280219856238}{https://doi.org/10.1177/0962280219856238}
}



\examples{
library(ALassoSurvIC)

### Getting the unpenalized NPMLE for interval censored data
data(ex_IC)
lowerIC <- ex_IC$lowerIC
upperIC <- ex_IC$upperIC
X <- ex_IC[, -c(1:2)]
system.time(result <- unpencoxIC(lowerIC, upperIC, X))

result           # main result
baseline(result) # obtaining the baseline cumulative hazard estimate
plot(result)     # plotting the estimated baseline cumulative hazard function by default
plot(result, what = "survival")  # plotting the estimated baseline survival function

## Use the parallel computing to reduce computational times
\donttest{library(parallel)
cl <- makeCluster(2L)  # making the cluster object 'cl' with two CPU cores
cl <- makeCluster(detectCores()) # run this code to use all CPU cores
system.time(result <- unpencoxIC(lowerIC, upperIC, X, cl = cl))}

\donttest{### Getting the unpenalized NPMLE for interval censored and left truncated data
## Try following codes with the 'ex_ICLT' data example
data(ex_ICLT)
lowerIC <- ex_ICLT$lowerIC
upperIC <- ex_ICLT$upperIC
trunc <- ex_ICLT$trunc
X <- ex_ICLT[, -c(1:3)]
result2 <- unpencoxIC(lowerIC, upperIC, X, trunc)
result2

baseline(result2)
plot(result2)
plot(result2, what = "survival")
}
}

\seealso{\code{alacoxIC}}
\keyword{functions}
