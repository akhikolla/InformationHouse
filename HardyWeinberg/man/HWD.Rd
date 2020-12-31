\name{HWD}
\alias{HWD}
\title{
  Compute disequilibrium statistic D 
}
\description{
Function \code{HWD} computes Weir's disequilibrium coefficient D.
}
\usage{
HWD(X)
}
\arguments{
  \item{X}{a vector of genotype counts (AA, AB, BB)}
}
\value{Returns the disequilibrium coefficient}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWf}} \code{\link{HWChisq}}}
\examples{
x <- c(MM=298,MN=489,NN=213)
D <- HWD(x)
cat("Disequilibrium coefficient: ",D,"\n")
}
\keyword{misc}
