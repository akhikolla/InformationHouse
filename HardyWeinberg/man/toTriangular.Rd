\name{toTriangular}
\alias{toTriangular}
\title{Convert a vector of genotype counts to triangular format}
\description{
Function \code{toTriangular} converts a named vector of genotype counts
into a triangular matrix format, with homozygotes on the diagonal and
heterozygotes below the diagonal.
}
\usage{
toTriangular(x)
}
\arguments{
  \item{x}{A vector of genotype counts}
}
\value{a matrix}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\examples{
x <- c(AA=20,AB=52,AC=34,BB=17,BC=51,CC=26)
print(x)
X <- toTriangular(x)
print(X)
}
\keyword{misc}

