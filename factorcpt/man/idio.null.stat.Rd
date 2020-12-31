\name{idio.null.stat}
\alias{idio.null.stat}

\title{
Threshold selection for the idiosyncratic components
}

\description{
Implements the stationary bootstrap procedure, which generates thresholds
that are used for multiple change-point detection in the idiosyncratic components.
It is internally called for by \code{\link{idio.seg}} and inherits all its arguments.
}

\usage{
idio.null.stat(mat, p, B, hat.vep, diag, scales, mby, tby, do.parallel)
}

\arguments{
  \item{mat}{
  a matrix containing the nodes of a binary tree grown for change-point analysis
  from the estimated idiosyncratic component
  }
  \item{p}{
  reciprocals of the average block size for the factors
  }
  \item{B}{
  the number of bootstrap samples
  }
  \item{hat.vep}{
  estimated idiosyncratic components
  }
  \item{diag}{
  see \code{idio.diag} in \code{\link{factor.seg.alg}}
  }
  \item{scales}{
  see \code{scales} in \code{\link{factor.seg.alg}}
  }
  \item{mby}{
  see \code{dmby} in \code{\link{func_dc_by}}
  }
  \item{tby}{
  see \code{dtby} in \code{\link{func_dc_by}}
  }
  \item{do.parallel}{
  see \code{do.parallel} in \code{\link{factor.seg.alg}}
  }
}

\value{
\code{B} values approximating the test statistics under the stationarity
over the segments of consideration in the binary tree (defined by the nodes contained in \code{mat}).
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint.
}

\author{
Haeran Cho
}

\keyword{internal}
