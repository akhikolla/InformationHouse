\name{ryan}
\alias{ryan}
\docType{data}
\title{Ryan data}
\description{
  This data set contains the data given in Table 9.2 by Ryan (2011, p. 323). The sample
  comprises 20 subgroups, each with 4 observations, on two quality characteristics \eqn{X_1}{X1} and
  \eqn{X_2}{X2}. According to Ryan (2011), the 10th and 20th subgroups are out-of-control.
}
\usage{data(ryan)}
\format{
  A 2x4x20 array.
}
\references{
  T. P. Ryan (2011), \emph{Statistical Methods for Quality Improvement},
  3rd ed., Wiley.
}
\examples{
data(ryan)
mphase1Plot(ryan)
}
\keyword{datasets}