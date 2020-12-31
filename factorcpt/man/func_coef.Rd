\name{func_coef}
\alias{func_coef}

\title{
Haar wavelet coefficients
}

\description{
Computes the Haar wavelet coefficients of input time series matrix at a given scale via filtering.
}

\usage{
func_coef(z, scale)
}

\arguments{
  \item{z}{
  input time series matrix
  }
  \item{scale}{
  a negative integer representing the wavelet scale
  }
}

\value{
A matrix containing Haar wavelet coefficients of each time series as its row
}

\author{
Haeran Cho
}

\keyword{internal}
