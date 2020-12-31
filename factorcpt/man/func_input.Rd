\name{func_input}
\alias{func_input}

\title{Haar wavelet-based transform}

\description{
Takes the Haar wavelet coefficients of multiple time series as an input and
applies the transforms proposed in Barigozzi, Cho & Fryzlewicz (2016), \code{g_j} and \code{h_j},
which turns \code{n}-dimensional time series to \code{n(n+1)/2}-dimensional panel data.
}

\usage{
func_input(coef, sgn)
}

\arguments{
  \item{coef}{
  Haar wavelet coefficients returned by \code{\link{func_coef}}
  }
  \item{sgn}{
  \code{n}-by-\code{n} matrix containing the sign of correlations between any pair of
  Haar wavelet coefficients in \code{coef}
  }
}

\value{
\code{n(n+1)/2}-dimensional panel of wavelet-transforms of input time series
}

\author{Haeran Cho}

\keyword{internal}
