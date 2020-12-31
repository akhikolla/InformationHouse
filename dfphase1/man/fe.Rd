\name{fe}
\alias{fe}
\docType{data}
\title{Ferric Oxide data}
\description{
  This data set contains 189 ferric-oxide individual measurement
  collected in an aluminum smelter.
}
\usage{data(fe)}
\format{
  A vector of length 189.
}
\references{
  M. D. Holland, D. M. Hawkins (2014) \dQuote{A Control Chart Based on a
    Nonparametric Multivariate Change-Point Model}, \emph{Journal of
    Quality Technology}, \bold{46}, 
  pp 63--77.
}
\examples{
data(fe)
phase1Plot(fe)
}
\keyword{datasets}