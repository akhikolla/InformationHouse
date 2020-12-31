\name{v_measure}
\alias{v_measure}

\title{Compute V-measure}
\description{Compute V-measure, a segmentation evaluation measure, which is based upon two criteria for clustering usefulness, homogeneity and completeness.}

\usage{
v_measure(sig, est, beta = 1)
}

\arguments{
\item{sig}{true signal; a numeric vector}
\item{est}{estimator; a numeric vector}
\item{beta}{parameter in definition of V-measure, see (Rosenberg and Hirschberg, 2007) for details}
}

\value{
A scalar takes value in [0, 1], with a larger value indicating higher accuracy. 
}

\references{
Rosenberg, A., and Hirschberg, J. (2007). V-measure: a conditional entropy-based external
cluster evaluation measures. \emph{Proc. Conf. Empirical Methods Natural Lang. Process.},
(June):410--420.
} 

\seealso{\code{\link{computeFdp}}, \code{\link{smuce}}, \code{\link{fdrseg}}, \code{\link{evalStepFun}}}

\examples{
# simulate data
u0 <- c(rep(1, 50), rep(5, 50))
Y  <- rnorm(100, u0)

# compute FDRSeg
uh <- fdrseg(Y)

plot(Y, pch = 20, col = "grey", xlab = "", ylab = "")
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "FDRSeg"), lty = c(1, 1), col = c("blue", "red"))

# compute V-measure
vm <- v_measure(u0, evalStepFun(uh))
print(vm)
}

\keyword{nonparametric}


