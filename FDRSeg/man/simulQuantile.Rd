\name{simulQuantile}
\alias{simulQuantile}

\title{Quantile simulations}
\description{Simulate the quantiles of multiscale statistics for SMUCE, FDRSeg, and D-FDRSeg under null hypothesis.}

\usage{
simulQuantile(alpha, n, r = round(50/min(alpha, 1-alpha)), 
              type = c("smuce","fdrseg","dfdrseg"), convKern, pos = .GlobalEnv)
}

\arguments{
\item{alpha}{a scalar with values in [0, 1]; the \code{alpha}-quantile of the null distribution of the multiscale statistic for SMUCE, FDRSeg, or D-FDRSeg via Monte Carlo simulation, see (Frick et al., 2014; Hotz et al., 2013; Li et al., 2015) for an explanation}
\item{n}{number of observations}
\item{r}{numer of Monte Carlo simulations}
\item{type}{
\preformatted{
"smuce"     simulate quantile for SMUCE
"fdrseg"    simulate quantiles for FDRSeg
"dfdrseg"   simulate quantiles for D-FDRSeg
}}
\item{convKern}{convolution kernel, only needed when \code{type} is {"dfdrseg"}}
\item{pos}{environment for saving the simulations for possible later usage}
}


\value{
A scalar value if \code{type} is chosen as {"smuce"}; a numeric vector of length \code{n} if \code{type} is chosen as {"fdrseg"} or {"dfdrseg"}.
}

\references{
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference. 
\emph{J. R. Statist. Soc. B, with discussion and rejoinder by the authors}, 76:495--580.

Hotz, T., Schuette, O. M., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013). Idealizing ion channel recordings by a jump segmentation multiresolution filter. \emph{IEEE Transactions on Nanobioscience}, 12(4):376--86. 

Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
} 

\seealso{\code{\link{smuce}}, \code{\link{fdrseg}}, \code{\link{dfdrseg}}}

\examples{
library(stepR)

# simulate quantiles for independent Gaussian noises
qs  <- simulQuantile(0.9, 100, type = "smuce")
qfs <- simulQuantile(0.9, 100, type = "fdrseg")
# plot result
yrng <- range(qs, qfs)
plot(qfs, pch = 20, ylim = yrng, xlab = "n", ylab = "")
abline(h = qs)

# simulate  quantiles for dependent Gaussian noises
convKern <- dfilter("bessel")$kern  # create digital filters
qdfs     <- simulQuantile(0.9, 100, type = "dfdrseg", convKern = convKern)
plot(qdfs, pch = 20, xlab = "n", ylab = "")
}

\keyword{nonparametric}


