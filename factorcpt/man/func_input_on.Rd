\name{func_input_on}
\alias{func_input_on}

\title{Haar wavelet-based diagonal transform}

\description{
Takes the Haar wavelet coefficients of multiple time series as an input and
applies the diagonal transforms proposed in Barigozzi, Cho & Fryzlewicz (2016), \code{g_j},
which turns \code{n}-dimensional time series to \code{n}-dimensional panel data.
}

\usage{
func_input_on(coef)
}

\arguments{
  \item{coef}{
  Haar wavelet coefficients returned by \code{\link{func_coef}}
  }
}

\value{
\code{n}-dimensional panel of wavelet-transforms of input time series
}

\author{
Haeran Cho
}

\keyword{internal}
