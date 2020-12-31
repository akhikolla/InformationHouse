\name{mac}
\alias{mac}
\title{
 Compute the minor allele count.
}
\description{
 \code{mac} computes the smallest allele count for a given vector of
 genotype counts.
}
\usage{
mac(X)
}
\arguments{
  \item{X}{a vector or matrix with genotype counts (AA, AB, BB)}
}
\value{
  a vector of the minor allele counts
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
   \code{\link{maf}}
}
\examples{
   X <- as.vector(rmultinom(1,100,c(0.5,0.4,0.1)))
   names(X) <- c("AA","AB","BB")
   print(X)
   print(mac(X))
}
\keyword{misc}