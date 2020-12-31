\name{HWExactMat}
\Rdversion{1.1}
\alias{HWExactMat}
\title{Matrix version of HWExact}
\description{
   HWExactMat executes a fast Exact test for HWE for each row in a matrix.
}
\usage{
HWExactMat(X, ...)
}
\arguments{
  \item{X}{A n times 3 matrix of genotypic counts (AA,AB,BB)}
  \item{\dots}{extra arguments that are passed on to \code{HWExact}}
}
\value{
  \item{pvalvec}{Vector with the p-values of each test}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWExact}} }
\examples{
X <- HWData(100,10)
colnames(X) <- c("MM","MN","NN")
Results <- HWExactMat(X)
Output <- cbind(X,Results$pvalvec)
print(Output)
}
\keyword{htest}
