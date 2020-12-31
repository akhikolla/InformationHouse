\name{HWChisqMat}
\Rdversion{1.1}
\alias{HWChisqMat}
\title{Matrix version of HWChisq}
\description{
   HWChisqMat executes the Chisquare test for HWE for each row in a matrix.
}
\usage{
HWChisqMat(X, ...)
}
\arguments{
  \item{X}{A n times 3 matrix of genotypic counts (AA,AB,BB)}
  \item{\dots}{extra arguments that are passed on to \code{HWChisq}}
}
\value{
  \item{pvalvec}{Vector with the p-values of each test}
  \item{chisqvec}{Vector with the chi-square statistics}
  \item{Dvec}{Vector with deviations from independence}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
X <- HWData(100,10)
colnames(X) <- c("MM","MN","NN")
Results <- HWChisqMat(X)
Output <- cbind(X,Results$chisqvec,Results$pvalvec)
print(Output)
}
\keyword{htest}
