\name{gravel}
\alias{gravel}
\docType{data}
\title{Gravel data}
\description{
  This data set contains 56 individual bivariate observations
  from a gravel-producing plant given by Holmes
  and Mergen (1993). There are two variables measuring
  the percentage of the particles (by weight) that are large or medium in size, respectively.
}
\usage{data(gravel)}
\format{
  A 2x56 matrix.
}
\references{
  D. S. Holmes, A. Mergen (1993) \dQuote{Improving the Performance of the \eqn{T^2}{T2} Control
    Chart}, \emph{Quality Engineering}, \bold{5}, pp. 619--625.
}
\examples{
data(gravel)
mphase1Plot(gravel)
}
\keyword{datasets}