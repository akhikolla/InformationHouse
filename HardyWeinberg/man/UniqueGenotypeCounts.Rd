\name{UniqueGenotypeCounts}
\Rdversion{1.1}
\alias{UniqueGenotypeCounts}
\title{Extract unique genotypic compositions from a matrix}
\description{
    Function \code{UniqueGenotypeCounts} creates a matrix containing
only the unique rows in the given matrix, together with their frequency
of occurrence
}
\usage{
UniqueGenotypeCounts(X, verbose = TRUE)
}
\arguments{
  \item{X}{A n by 3 matrix with genotypic counts (AA,AB,BB)}
  \item{verbose}{If TRUE then print some statistics}
}
\value{
  A matrix with 4 columns, AA, AB, BB, and frequency of occurrence
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{GenerateSamples}} }
\examples{
set.seed(123)
X <- HWData(n=100,nm=100)
print(nrow(X))
Y <- UniqueGenotypeCounts(X)
print(nrow(Y))
print(sum(Y$w))
}
\keyword{misc}

