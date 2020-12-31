\name{plot.clogitboost}
\alias{plot.clogitboost}

\title{
Plotting after fitting a boosting conditional logit model
}

\description{
  \code{plot} methods for the \code{clogitboost} objects, which produce marginal plots of the covariate effects.
}
  
\usage{
\method{plot}{clogitboost}(x, d, grid = NULL, \dots)
}

\arguments{
    \item{x}{output object from the \code{\link{clogitboost}} function.
        }
    \item{d}{integer indicating which covariate is used.}
\item{grid}{grid of values for plotting. If it is not specified, the minimal and maximal elements of the covariate are used as the two endpoints of the grid.}
    \item{\dots}{other options for plotting.}
}
\author{
Haolun Shi \email{shl2003@connect.hku.hk}

Guosheng Yin \email{gyin@hku.hk}
}

\seealso{
\code{\link{clogitboost}}
}
\examples{
data(travel)
train <- 1:504
y <- travel$MODE[train]
x <- travel[train, 3:6]
strata <- travel$Group[train]
fit <- clogitboost(y = y, x = x, strata = strata, iter = 10, rho = 0.05)
plot(fit, d = 1, xlab = "x", ylab = "f(x)", main = "TTIME", type = "l")

}
