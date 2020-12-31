\name{alleles}
\alias{alleles}
\title{
Extract alleles 
}
\description{
Function \code{alleles} extracts the names of the alleles from a named
genotype vector.
}
\usage{
alleles(x, fromlabels = TRUE)
}
\arguments{
  \item{x}{A named or unnamed genotype vector (e.g. c(AA=10,AB=20,BB=5))}
  \item{fromlabels}{extract genotypes from the labels of the vector
    elements, or from the vector elements themselves.} 
}
\value{A character vector with the alleles}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{\code{\link{n.alleles}}}
\examples{
   x <- c(AA=10,AG=10,GG=10,AT=5)
   als.x <- alleles(x)
}
\keyword{misc}
