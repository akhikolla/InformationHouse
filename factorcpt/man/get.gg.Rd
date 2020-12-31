\name{get.gg}
\alias{get.gg}

\title{
Data-driven selection of the average block size
}

\description{
Computes quantities required for data-driven selection of the average block size for stationary bootstrap
}

\usage{
get.gg(z, M = NULL, C = 2, max.K = 5)
}

\arguments{
  \item{z}{
  univariate time series
  }
  \item{M}{
  bandwidth, by default, \code{M = NULL} to be automatically selected as in Politis and White (2004)
  }
  \item{C}{
  \code{C = 2} is the default value chosen for automatic bandwidth selection in Politis and White (2004)
  }
  \item{max.K}{
  \code{max.K = 5} is the default value chosen for automatic bandwidth selection in Politis and White (2004)
  }
}

\value{
Estimates for the two quantities required for average block size selection,
see the Appendix of Barigozzi, Cho & Fryzlewicz (2016).
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint. \cr

D. N. Politis and H. White (2004) Automatic block-length selection for the dependent bootstrap. Econometric Reviews. 23: 53-70.
}

\author{Haeran Cho}

\keyword{internal}
