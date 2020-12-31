\name{HWAlrPlot}
\alias{HWAlrPlot}
\title{Plot genetic markers in additive log-ratio coordinates}
\description{
  \code{HWAlrPlot} creates a scatter plot of the log-ratio coordinates of
  bi-allelic genetic markers. Hardy-Weinberg equilibrium is indicated by
  a straigh line in the plot.
}
\usage{
HWAlrPlot(X, zeroadj = 0.5)
}
\arguments{
  \item{X}{A matrix of genotype counts (columns AA, AB, BB)
}
  \item{zeroadj}{Zero-adjustment parameter. Zero counts in the count
    matrix are substituted by \code{zeroadj} which is 0.5 by default.
}
}
\value{
   NULL
 }
\references{Graffelman, J. and Egozcue, J. J. (2011) Hardy-Weinberg
  equilibrium: a non-parametric compositional approach. In: Vera
  Pawlowsky-Glahn and Antonella Buccianti (eds.) Compositional Data
  Analysis: Theory and Applications, John Wiley \& Sons, Ltd, pp. 207-215}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{HWClrPlot},\link{HWIlrPlot}}
}
\examples{
   X <- HWClo(HWData(100,100))
   HWAlrPlot(X)
}
\keyword{aplot}