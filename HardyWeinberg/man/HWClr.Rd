\name{HWClr}
\alias{HWClr}
\title{
  Compute the centred log-ratio transformation
}
\description{
  \code{HWClr} computes the centred log-ratio transformation for
  genotype counts of bi-allelic genetic markers.
}
\usage{
HWClr(X, zeroadj = 0.5)
}
\arguments{
  \item{X}{A matrix of genotype counts (columns AA, AB and BB)}
  \item{zeroadj}{A zero adjustment parameter (0.5 by default)}
}
\value{
  A matrix or vector of log-ratio coordinates
}
\references{Graffelman, J. and Egozcue, J. J. (2011) Hardy-Weinberg
  equilibrium: a non-parametric compositional approach. In: Vera
  Pawlowsky-Glahn and Antonella Buccianti (eds.) Compositional Data
  Analysis: Theory and Applications, John Wiley \& Sons, Ltd, pp. 207-215}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{HWAlr},\link{HWIlr}}
}
\examples{
   X <- HWData(100,100)
   Y <- HWClr(X)
}
\keyword{misc}
