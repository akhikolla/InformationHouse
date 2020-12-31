\name{HWf}
\alias{HWf}
\title{
 Computation of inbreeding coefficient
}
\description{
  \code{HWf} computes the inbreeding coefficient for sample
  of genotype counts, or a matrix of genotype counts.
}
\usage{
HWf(X)
}
\arguments{
  \item{X}{a vector or matrix of genotype counts (AA, AB, BB)}
}
\details{
   For monomorphic markers a warning is issued, and the estimate for the
   inbreeding coefficient is NaN. 
}
\value{
  Returns a single inbreeding coefficient (intraclass correlation coefficient), if \code{X} is a single sample, or a vector of inbreeding coefficients, if \code{X} is a matrix with genotype counts.
}
\references{
Crow, J. F. and Kimura, M. (1970) An introduction to population genetics theory.
Harper \& Row, publishers, New York
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
#
# A single sample
#
x <- c(MM=298,MN=489,NN=213)
fhat <- HWf(x)
cat("Inbreeding coefficient: ",fhat,"\n")

#
# Multiple samples
#
X <- HWData(nm=100,n=1000)
fhat <- HWf(X)
}
\keyword{misc}




