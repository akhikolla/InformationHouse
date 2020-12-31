\name{n.alleles}
\alias{n.alleles}
\title{
Number of alleles 
}
\description{
Function \code{n.alleles} determines the number of alleles in a named
genotype vector.
}
\usage{
n.alleles(x, ...)
}
\arguments{
  \item{x}{A named genotype vector (e.g. c(AA=10,AB=20,BB=5))}
  \item{\dots}{extra arguments that are passed on to \code{alleles}}
}
\value{integer}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{\code{\link{alleles}}}
\examples{
   x <- c(AA=25,AB=50,BB=25)
   k <- n.alleles(x)
   print(k)
}
\keyword{misc}

