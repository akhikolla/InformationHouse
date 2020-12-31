\name{HWAlr}
\alias{HWAlr}
\title{
  Compute additive log-ratio transformation
}
\description{
  \code{HWAlr} computes the additive log-ratio transformation for
  genotype counts of bi-allelic genetic markers.
}
\usage{
HWAlr(X, zeroadj = 0.5, denominator = 2)
}
\arguments{
  \item{X}{A matrix of genotype counts (columns AA, AB and BB)}
  \item{zeroadj}{A zero adjustment parameter (0.5 by default)}
  \item{denominator}{The genotype count put in the denominator of the
    log-ratio (1=AA, 2=AB, 3=BB)}
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
\code{\link{HWClr},\link{HWIlr}}
}
\examples{
   X <- HWData(100,100)
   Y <- HWAlr(X)
}
\keyword{misc}

