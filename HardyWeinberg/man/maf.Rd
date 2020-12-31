\name{maf}
\alias{maf}
\title{Function to compute minor allele frequencies}
\description{
  Function \code{maf} computes the minor allele frequency for
  a matrix or vector of compositions.
}
\usage{
maf(x)
}
\arguments{
  \item{x}{a vector or matrix of genotypic compositions}
}
\value{
  a vector of minor allele frequencies.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
   X <- as.vector(rmultinom(1,100,c(0.5,0.4,0.1)))
   X <- X/sum(X)
   print(X)
   print(maf(X))
}
\keyword{misc}

