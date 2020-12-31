\name{af}
\alias{af}
\title{Function to compute allele frequencies}
\description{
  Function \code{af} computes the allele frequencies for
  a matrix or a vector containing genotypic compositions.
}
\usage{
af(x)
}
\arguments{
  \item{x}{a vector or matrix with compositions}
}
\value{
  a vector with allele frequencies
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{maf}}}
\examples{
   X <- as.vector(rmultinom(1,100,c(0.5,0.4,0.1)))
   X <- X/sum(X)
   print(X)
   print(af(X))
}
\keyword{misc}

