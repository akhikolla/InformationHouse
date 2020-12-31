\name{smuce}
\alias{smuce}

\title{Piecewise constant regression with SMUCE}
\description{Compute the SMUCE estimator for one-dimensional data with i.i.d. Gaussian noises.}

\usage{
smuce(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)), sd = stepR::sdrobnorm(Y))
}

\arguments{
\item{Y}{a numeric vector containing the noisy data}
\item{q}{threshold value; a scalar number}
\item{alpha}{significance level; if \code{q} is missing, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see (Frick et al., 2014) for an explanation}
\item{r}{numer of Monte Carlo simulations}
\item{sd}{standard deviation of noises}
}

\value{
A {list} with components
\item{value}{function values on each segment of the estimator}
\item{left}{indices of leftmost points within each segment of the estimator}
\item{n}{number of samples}
}

\note{
This is an efficient implementation of function \code{smuceR} in R package stepR (CRAN) for data with i.i.d. Gaussian noises. The detailed algorithm is described in (Seiling, 2013).
}

\references{
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference. 
\emph{J. R. Statist. Soc. B, with discussion and rejoinder by the authors}, 76:495--580.

Seiling, H. (2013). Statistical Multiscale Segmentation: Inference, Algorithms and Applications. 
PhD thesis, University of Goettingen, Germany.
} 

\seealso{\code{\link{fdrseg}}, \code{\link{dfdrseg}}, \code{\link{simulQuantile}}, \code{\link{sdrobnorm}}, \code{\link{evalStepFun}}, \code{\link{computeFdp}}, \code{\link{v_measure}}}

\examples{
# simulate data
set.seed(2)
u0 <- c(rep(1, 50), rep(5, 50))
Y  <- rnorm(100, u0)

# compute the estimate (q is automatically simulated)
uh <- smuce(Y)

# plot result
plot(Y, pch = 20, col = "grey", ylab = "", main = expression(alpha*" = 0.1"))
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "SMUCE"), lty = c(1, 1), col = c("blue", "red"))


# other choice of alpha
uh <- smuce(Y, alpha = 0.05)
# plot result
plot(Y, pch = 20, col = "grey", ylab = "", main = expression(alpha*" = 0.05"))
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "SMUCE"), lty = c(1, 1), col = c("blue", "red"))

\dontrun{
# alternatively simulate quantiles first
alpha <- 0.1
q     <- simulQuantile(1 - alpha, 100, type = "smuce")

# then compute the estimate
uh <- smuce(Y, q)}}

\keyword{nonparametric}
