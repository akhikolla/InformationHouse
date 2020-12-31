\name{dlevene}
\alias{dlevene}
\title{
Calculate Levene's exact density for k alleles
}
\description{
Function \code{dlevene} calculates Levene's exact density for a diploid
system with k alleles.
}
\usage{
dlevene(N)
}
\arguments{
  \item{N}{A lower triangular matrix with genotype counts}
}
\details{
  The supplied matrix of genotype counts should be triangular, with the
  homozygote counts on the diagonal, and all heterozygote counts below
  the diagonal.
}
\value{a single real number}
\references{
Levene, H. (1949) On a matching problem arising in genetics. Annals of
Mathematical Statistics, 20,  pp. 91-94.
}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{
\code{\link{HWExact}}
}
\examples{
x <- c(AA=12,AB=19,AC=13,BB=7,BC=5,CC=0)
x <- toTriangular(x)
prob <- dlevene(x)
print(prob)
}
\keyword{distribution}

