\name{func_dc_by}
\alias{func_dc_by}

\title{
Double CUSUM statistics
}

\description{
Aggregates high-dimensional CUSUM statistic series via Double CUSUM approach,
with additional parameters introduced to speed up the computation.
}

\usage{
func_dc_by(z, dmby, dtby)
}

\arguments{
  \item{z}{
  input data matrix
  }
  \item{dmby}{
  the function computes every \code{dmby}th Double CUSUM statistic cross-sectionally
  }
  \item{dtby}{
  the function computes every \code{dtby}th Double CUSUM statistic in time
  }
}

\value{
\item{cs}{a matrix of CUSUM series}
\item{acs}{a matrix of absolute values of CUSUM series}
\item{res}{pointwise maximum of Double CUSUM statistics}
\item{mat}{a matrix of Double CUSUM statistics}
}

\references{
H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. Electronic Journal of Statistics. 10: 2000-2038.}

\author{
Haeran Cho
}

\keyword{internal}
