\name{HWIlrPlot}
\alias{HWIlrPlot}
\title{
Plot bi-allelic genetic markers in isometric log ratio coordinates
}
\description{
\code{HWIlrPlot} makes a scatter plot of the isometric log ratio
coordinates for bia-llelic markers.
}
\usage{
HWIlrPlot(X, zeroadj = 0.5, ...)
}
\arguments{
  \item{X}{Matrix of genotype counts, one marker per row, AA, AB and BB
    in three columns
}
  \item{zeroadj}{Adjustment for zero values (0.5 by default)
}
  \item{\dots}{Additional arguments for function \code{plot}
}
}
\value{
A matrix of log ratio coordinates.
}
\references{
  Graffelman, J. and Egozcue, J. J. (2011) Hardy-Weinberg equilibrium: a
  nonparametric compositional approach. 
  In Pawlowsky-Glahn, V. and Buccianti A., editors, Compositional Data Analysis: Theory and Applications, pages 208-217, John Wiley & Sons, Ltd. 
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{HWAlrPlot},\link{HWClrPlot}}}
\examples{
X <- HWClo(HWData(100,100))
HWIlrPlot(X)
}
\keyword{aplot}
