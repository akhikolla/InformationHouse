\name{genlabels}
\alias{genlabels}
\title{
Label genotype counts of a vector or matrix
}
\description{
Function \code{genlabels} sets the names of a vector or matrix of
genotype counts.
}
\usage{
genlabels(X)
}
\arguments{
  \item{X}{a 3 (or 5) element vector with genotype counts, a matrix of
    genotype counts (3 or 5 columns)}
}
\value{A vector or a matrix}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{HWChisq}}}
\examples{
x <- c(25,50,25)
x <- genlabels(x)
}
\keyword{~misc}
