\name{plot.voronoi.polygons}
\alias{plot.voronoi.polygons}
\title{plots an voronoi.polygons object}
\description{
  plots an \code{voronoi.polygons} object
}
\usage{
\method{plot}{voronoi.polygons}(x, which, color=TRUE, isometric=TRUE, ...)
}
\arguments{
  \item{x}{ object of class \code{voronoi.polygons}  }
  \item{which}{ index vector selecting which polygons to plot }
  \item{color}{ logical, determines if plot should be colored, default:
    \code{TRUE} }
  \item{isometric}{generate an isometric plot (default \code{TRUE})}
  \item{\dots}{ additional plot arguments }
}
\author{ A. Gebhardt}
\seealso{ \code{\link{voronoi.polygons}}}
\examples{
data(franke)
fd3 <- franke$ds3
fd3.vm <- voronoi.mosaic(fd3$x,fd3$y)
fd3.vp <- voronoi.polygons(fd3.vm)
plot(fd3.vp)
plot(fd3.vp,which=c(3,4,6,10))
}
\keyword{ spatial }% at least one, from doc/KEYWORDS
