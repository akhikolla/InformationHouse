\name{HWLindley}
\alias{HWLindley}
\title{
  Calculate a posteriori density for Lindley's alpha
}
\description{
  Function \code{HWLindley} calculates the posterior density for disequilibrium measure alpha, as defined by Lindley (1988).
}
\usage{
HWLindley(alphaseq = seq(-3, 3, by = 0.01), x)
}
\arguments{
  \item{alphaseq}{a single value or a sequence of values for alpha}
  \item{x}{the genotype count vector in format (AA,AB,BB)}
}
\details{
  Numerical integration is used to compute the density.
}
\value{
  a vector with values of the density for each value in \code{alphaseq}
}
\references{
  Lindley, D.V. (1988) Statistical Inference Concerning Hardy-Weinberg Equilibrium. In: 
  Bernardo, J.M., DeGroot, M.H., Lindley, D.V. and Smith, A.F.M. Bayesian Statistics, 3,
  pp. 307-326. Oxford University Press.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
   \code{\link{HWPosterior}}
}
\examples{
x <- c(MM=298,MN=489,NN=213)
post.dens <- HWLindley(seq(-2,2,by=0.01),x)
\dontrun{
  plot(seq(-2,2,by=0.01),post.dens,type="l")
}
}
\keyword{distribution}

