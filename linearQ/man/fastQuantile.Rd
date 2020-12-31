\name{fastQuantile}
\alias{fastQuantile}
%\docType{package}
\title{
linear Algorithm for Quantile Simulation
}
\description{
This is a linear algorithm for quantile simulation under null hypothesis in multiscale change-point segmentation.
}
\usage{
fastQuantile(alpha, n, r=round(50/min(alpha, 1-alpha)), 
              mType=c("norm-pen","pois"), seed = 123, ...)	
}

\arguments{
\item{alpha}{a scalar with values in [0, 1]; the \code{alpha}-quantile of the null distribution of the multiscale statistic via Monte Carlo simulation}
\item{n}{number of observations}
\item{r}{number of Monte Carlo simulations}
\item{mType}{"norm-pen" simulates the multiscale statistic from Normal regression model, "pois" simulates the multiscale statistic from Poission regression model. }
\item{seed}{data seed}
\item{...}{further arguments passed to penalty function}
}

\value{
 A scalar quantile value q.
}

\references{
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference.
\emph{J. R. Statist. Soc. B, with discussion and rejoinder by the authors}, 76:495--580.

Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
}

\examples{
# simulate quantiles for multiscale statistics from Normal regression model
    seed = 123
    q  <- fastQuantile(0.9, 500, 100, mType = "norm-pen")
}
