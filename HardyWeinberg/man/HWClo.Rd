\name{HWClo}
\alias{HWClo}
\title{
Convert genotype counts to compositions
}
\description{
Function \code{HWClo} divides each row of a matrix by its total, and so
produces matrix of compositions.
}
\usage{
HWClo(X)
}
\arguments{
  \item{X}{A matrix of (genotype) counts}
}
\value{
A matrix
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
  \code{\link{HWAlr},\link{HWClr},\link{HWIlr}}
}
\examples{
  X <- HWData(2,100)
  Y <- HWClo(X)
}
\keyword{misc}
