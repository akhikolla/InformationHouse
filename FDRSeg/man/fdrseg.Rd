\name{fdrseg}
\alias{fdrseg}

\title{Piecewise constant regression with FDRSeg}
\description{Compute the FDRSeg estimator for one-dimensional data with i.i.d. Gaussian noises.}

\usage{
fdrseg(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)), sd = stepR::sdrobnorm(Y))
}

\arguments{
\item{Y}{a numeric vector containing the noisy data}
\item{q}{threshold value; a numeric vector of the same length as the data}
\item{alpha}{significance level; if \code{q} is missing, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see (Li et al., 2015) for an explanation}
\item{r}{numer of Monte Carlo simulations}
\item{sd}{standard deviation of noises}
}

\value{
A {list} with components
\item{value}{function values on each segment of the estimator}
\item{left}{indices of leftmost points within each segment of the estimator}
\item{n}{number of samples}
}

\references{
Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
} 

\seealso{\code{\link{smuce}}, \code{\link{dfdrseg}}, \code{\link{simulQuantile}}, \code{\link{sdrobnorm}}, \code{\link{evalStepFun}}, \code{\link{computeFdp}}, \code{\link{v_measure}}}

\examples{
# simulate data
set.seed(123)
u0 <- c(rep(1, 50), rep(5, 50))
Y  <- rnorm(100, u0)

# compute the estimate (q is automatically simulated)
# it might take a while due to simulating quantiles and will
# be faster for later calls on signals of the same length
uh <- fdrseg(Y)

# plot result
plot(Y, pch = 20, col = "grey", ylab = "", main = expression(alpha*" = 0.1"))
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "FDRSeg"), lty = c(1, 1), col = c("blue", "red"))


# other choice of alpha
uh <- fdrseg(Y, alpha = 0.05)
# plot result
plot(Y, pch = 20, col = "grey", ylab = "", main = expression(alpha*" = 0.05"))
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "FDRSeg"), lty = c(1, 1), col = c("blue", "red"))

\dontrun{
# alternatively simulate quantiles first
alpha <- 0.1
q     <- simulQuantile(1 - alpha, 100, type = "fdrseg")

# then compute the estimate
uh <- fdrseg(Y, q)}
}

\keyword{nonparametric}

