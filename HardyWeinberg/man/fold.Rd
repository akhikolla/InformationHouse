\name{fold}
\alias{fold}
\title{
  Fold a square matrix
}
\description{
  The function \code{fold} sums corresponding below and above diagonal elements of a square matrix to form a triangular matrix. 
}
\usage{
fold(X, lower = TRUE)
}
\arguments{
  \item{X}{a square matrix}
  \item{lower}{logical. If \code{lower=TRUE} a lower triangular matrix is formed, if not an upper triangular matrix.
}
}
\details{
  Useful for constructing triangular matrices of genotype counts
}
\value{
  A matrix
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
\code{\link{lower.tri},\link{upper.tri}}
}
\examples{
allelenames <- paste("A",11:13,sep="")
allele1 <- factor(c("A11","A11","A12","A12","A13","A12"),levels=allelenames)
allele2 <- factor(c("A11","A12","A12","A13","A13","A11"),levels=allelenames)
GC  <- table(allele1,allele2)
GC  <- as.matrix(unclass(GC))
GCf <- fold(GC)
}
\keyword{array}
\keyword{manip}

