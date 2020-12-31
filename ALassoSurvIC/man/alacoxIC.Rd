\name{alacoxIC}
\alias{alacoxIC}
\alias{alacoxIC.default}
\alias{print.alacoxIC}
\title{Performing variable selection with an adaptive lasso penalty for interval censored and possibly left truncated data}
\usage{
\method{alacoxIC}{default}(lowerIC, upperIC, X, trunc, theta,
  normalize.X = TRUE, cl = NULL, max.theta = 1000, tol = 0.001,
  niter = 1e+05, string.cen = Inf, string.missing = NA, ...)
}
\arguments{
\item{...}{for S4 method only.}

\item{lowerIC}{A numeric vector for the lower limit of the censoring interval.}

\item{upperIC}{A numeric vector for the upper limit of the censoring interval.}

\item{X}{A numeric matrix for the covariates that will be used for variable selection.}

\item{trunc}{A numeric vector for left truncated times. If supplied, the function performs the variable selection for interval censored and left truncated data. If \code{trunc} is missing, the data will be considered as interval censored data.}

\item{theta}{A numeric value for the thresholding parameter. If \code{theta} is missing, the function automatically determines the thresholding parameter using a grid search algorithm, based on the Bayesian information criterion (BIC). See details below.}

\item{normalize.X}{A logical value: if \code{normalize.X} = \code{TRUE}, the covariate matrix \code{X} will be normalized before fitting models. Default is \code{TRUE}.}

\item{cl}{A cluster object created by \code{\link[parallel]{makeCluster}} in the \code{parallel} package. If \code{NULL}, no parallel computing is used by default. See details below.}

\item{max.theta}{A numeric value for the maximum value that a thresholding parameter can take when searching the optimal one. The algorithm will look up an optimal tunning parameter below \code{max.theta}. See details below.}

\item{tol}{A numeric value for the absolute iteration convergence tolerance.}

\item{niter}{A numeric value for the maximum number of iterations.}

\item{string.cen}{A string indicating right censoring for \code{upperIC}. Default is \code{Inf}.}

\item{string.missing}{A string indicating missing value. Default is \code{NA}.}

}
\description{
The \code{alacoxIC} function performs variable selection with an adaptive lasso penalty for interval censored and possibly left truncated data. It performs penalized nonparametric maximum likelihood estimation through a penalized EM algorithm by following Li et al. (2019). The function searches the optimal thresholding parameter automatically, based on BIC. The variable selection approach, implemented by the \code{alacoxIC} function, is proven to enjoy the desirable oracle property introduced by Fan & Li (2001). The full details are available in Li et al. (2019).
}

\details{
The grid search algorithm is used to find the optimal thresholding parameter using a grid search algorithm, based on BIC. Specifically, the \code{alacoxIC} function first searches the smallest integer thresholding parameter which all coefficient estimates are zero beween \eqn{1} and \code{max.theta} and then creates one hundred grid points by following the rule of Simon et al. (2011, Section 2.3). The one minimizing BIC among the one hundred candidates is chosen as the optimal thresholding parameter in the adaptive lasso estimation.

The cluster object, created by \code{\link[parallel]{makeCluster}} in the \code{parallel} package, can be supplied with the \code{cl} argument to reduce computation time via parallel computing. The parallel computing will be used when searching the optimal thresholding parameter and calculating the hessian matrix of the log profile likelihood. How to use the parallel computing is illustrated in one of the examples given below.

Use the \code{baseline} function and the \code{plot} function to extract and plot the estimate of the baseline cumulative hazard function, respectively, from the object returned by the \code{alacoxIC}. The \code{plot} function also provides the plot of the estimated baseline survival function. See the usages in the examples given below.
}

\references{
Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. Journal of the American statistical Association, 96(456), 1348-1360

Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2011). Regularization paths for Cox’s proportional hazards model via coordinate descent. Journal of statistical software, 39(5), 1.

Li, C., Pak, D., & Todem, D. (2019). Adaptive lasso for the Cox regression with interval censored and possibly left truncated data. Statistical methods in medical research. \href{https://doi.org/10.1177/0962280219856238}{https://doi.org/10.1177/0962280219856238}
}

\examples{
library(ALassoSurvIC)

### Variable selection for interval censored data
data(ex_IC) # the 'ex_IC' data having 100 subjects and 6 covariates
lowerIC <- ex_IC$lowerIC
upperIC <- ex_IC$upperIC
X <- ex_IC[, -c(1:2)]

## Performing the variable selection algorithm using a single core
\donttest{system.time(result <- alacoxIC(lowerIC, upperIC, X))}

## Use parallel computing to reduce the computation time
library(parallel)
cl <- makeCluster(2L)  # making the cluster object 'cl' with two CPU cores
\donttest{cl <- makeCluster(detectCores()) # run this code to use all CPU cores}

system.time(result <- alacoxIC(lowerIC, upperIC, X, cl = cl))

result           # main result
baseline(result) # obtaining the baseline cumulative hazard estimate
plot(result)     # plotting the baseline estimated cumulative hazard function by default
plot(result, what = "survival")  # plotting the estimated baseline survival function

### Variable selection for interval censored and left truncated data
## Try following codes with the 'ex_ICLT' data example
\donttest{data(ex_ICLT) # the 'ex_ICLT' data having 100 subjects and 6 covariates
lowerIC <- ex_ICLT$lowerIC
upperIC <- ex_ICLT$upperIC
trunc <- ex_ICLT$trunc
X <- ex_ICLT[, -c(1:3)]
result2 <- alacoxIC(lowerIC, upperIC, X, trunc)
result2

baseline(result2)
plot(result2)
plot(result2, what = "survival")
}
}

\seealso{\code{unpencoxIC}}
\keyword{functions}
