\name{computeFdp}
\alias{computeFdp}

\title{Compute false discovery proportion (FDP)}
\description{Compute false discovery proportion for estimated change-points, see (Li et al., 2015) for a detailed explanation.}

\usage{
computeFdp(u, eJ)
}

\arguments{
\item{u}{true signal; a numeric vector}
\item{eJ}{estimated change-points; a numeric vector}
}

\value{
A scalar takes value in [0, 1]. 
}

\references{
Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
} 

\seealso{\code{\link{fdrseg}}, \code{\link{v_measure}}}

\examples{
# simulate data
set.seed(2)
u0 <- c(rep(1, 50), rep(5, 50))
Y  <- rnorm(100, u0)

# compute FDRSeg
uh <- fdrseg(Y)

plot(Y, pch = 20, col = "grey", xlab = "", ylab = "")
lines(u0, type = "s", col = "blue")
lines(evalStepFun(uh), type = "s", col = "red")
legend("topleft", c("Truth", "FDRSeg"), lty = c(1, 1), col = c("blue", "red"))

# compute false discovery proportion
fdp <- computeFdp(u0, uh$left)
cat("False discovery propostion is ", fdp, "\n")
}

\keyword{nonparametric}
