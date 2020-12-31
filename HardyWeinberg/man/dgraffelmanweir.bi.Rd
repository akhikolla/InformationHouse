\name{dgraffelmanweir.bi}
\alias{dgraffelmanweir.bi}
\title{
Calculate Graffelman-Weir exact density for bi-allelic X-chromosomal variant
}
\description{
Function \code{dgraffelmanweir} calculate the probability P(NAB=nab and
MA=ma|NA=na) for a bi-allelic X-chromosomal variant.
}
\usage{dgraffelmanweir.bi(x, y)}
\arguments{
  \item{x}{vector with male genotype counts (A,B)}
  \item{y}{vector with female genotype counts (AA,AB,BB)}
}
\value{a single real number}
\references{
  Graffelman, J. and Weir, B.S. (2016) Testing for Hardy-Weinberg
  equilibrium at bi-allelic genetic markers on the X chromosome.
  \href{http://dx.doi.org/10.1038/hdy.2016.20}{Heredity} 116(6) pp. 558--568.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu}}
\seealso{\code{\link{HWExact}}, \code{\link{HWExactStats}}}
\examples{
males <- c(A=392, B=212)
females <- c(AA=275, AB=296, BB=80)
prob <- dgraffelmanweir.bi(males,females)
print(prob)
}
\keyword{distribution}

