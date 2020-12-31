\name{ThetatoF}
\alias{ThetatoF}
\title{
  Convert theta to an inbreeding coefficient
}
\description{
  Function \code{ThetatoF} converts disequilibrium measure theta to an
  inbreeding coefficient.
}
\usage{
ThetatoF(p, theta = 4)
}
\arguments{
  \item{p}{the allele frequency}
  \item{theta}{the disequilibrium parameter}
}
\value{
  a real number
}
\references{Rohlfs, R.V. and Weir, B.S. (2008) Distributions of
  {H}ardy-{W}einberg Equilibrium Test Statistics, Genetics, 180(3) pp. 1609-1616.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
   \code{\link{HWf}}
}
\examples{
f <- ThetatoF(0.5,4)
}
\keyword{misc}
