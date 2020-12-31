\name{dfdrseg}
\alias{dfdrseg}
\title{
Piecewise constant regression with D-FDRSeg
}
\description{
Compute the D-FDRSeg estimator for one-dimensional data with dependent Gaussian noises, especially for ion channel recordings, see (Hotz et al., 2013; Li et al., 2015) for further details.
}
\usage{
dfdrseg(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)), convKern, 
        sd = stepR::sdrobnorm(Y, lag=length(convKern)+1))
}

\arguments{
  \item{Y}{a numeric vector containing the noisy data}
  \item{q}{threshold value; a numeric vector of the same length as the data}
  \item{alpha}{significance level; if \code{q} is missing, \code{q} is chosen as the (1-\code{alpha})-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation, see (Li et al., 2015) for an explanation}
  \item{r}{numer of Monte Carlo simulations}
  \item{convKern}{kernel of the low-pass filter, see (Li et al., 2015)}
  \item{sd}{standard deviation of noises}
}

\value{
A {list} with components
\item{value}{function values on each segment of the estimator}
\item{left}{indices of leftmost points within each segment of the estimator}
\item{n}{number of samples}
}
\references{
Hotz, T., Schuette, O. M., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013). Idealizing ion channel recordings by a jump segmentation multiresolution filter. \emph{IEEE Transactions on Nanobioscience}, 12(4), 376-86. 

Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
}

\seealso{
\code{\link{smuce}}, \code{\link{dfdrseg}}, \code{\link{jsmurf}}, \code{\link{simulQuantile}}, \code{\link{sdrobnorm}}, \code{\link{contMC}}, \code{\link{dfilter}}, \code{\link{evalStepFun}}
}
\examples{
library(stepR)

# simulate data (a continuous time Markov chain)
ts        <- 0.1  # sampling time
SNR       <- 3    # signal-to-noise ratio
sampling  <- 1e4  # sampling rate 10 kHz
over      <- 10   # tenfold oversampling
cutoff    <- 1e3  # 1 kHz 4-pole Bessel-filter, adjusted for oversampling
simdf     <- dfilter("bessel", list(pole=4, cutoff=cutoff/sampling/over))
transRate <- 50
rates     <- rbind(c(0, transRate), c(transRate, 0))
set.seed(123)
sim <- contMC(ts*sampling, c(0,SNR), rates, sampling = sampling, family = "gaussKern",  
              param = list(df=simdf, over=over, sd=1))
Y   <- sim$data$y
x   <- sim$data$x

# D-FDRseg 
library(stepR)
convKern <- dfilter("bessel", list(pole=4, cutoff=cutoff/sampling))$kern
uh       <- dfdrseg(Y, convKern = convKern, r = 10) # r could be much larger

# plot results
plot(x, Y, pch = 20, col = "grey", xlab="", ylab = "", main = "Simulate Ion Channel Data")
lines(sim$discr, col = "blue")
lines(x, evalStepFun(uh), col = "red")
legend("topleft", c("Truth", "D-FDRSeg"), lty = c(1, 1), col = c("blue", "red"))

\dontrun{
# alternatively simulate quantiles first
alpha <- 0.1
q     <- simulQuantile(1 - alpha, ts*sampling, type = "dfdrseg", convKern = convKern)

# then compute the estimate
uh <- dfdrseg(Y, q, convKern = convKern)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
