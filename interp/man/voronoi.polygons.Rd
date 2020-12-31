\name{voronoi.polygons}
\alias{voronoi.polygons}
\title{ extract polygons from a voronoi mosaic }
\description{
  This functions extracts polygons from a \code{voronoi.mosaic} object.
}
\usage{
voronoi.polygons(voronoi.obj)
}
\arguments{
  \item{voronoi.obj}{ object of class \code{voronoi.mosaic} }
}
\value{
  Returns an object of class \code{voronoi.polygons} with unamed list
  elements for each polygon. These list
  elements are matrices with columns \code{x} and \code{y}.
  Unbounded polygons along the border are represented by \code{NULL}
  instead of a matrix.
}
\author{ Denis White }
\seealso{ \code{\link{plot.voronoi.polygons}},\code{\link{voronoi.mosaic}}}

\examples{
data(franke)
fd3 <- franke$ds3
fd3.vm <- voronoi.mosaic(fd3$x,fd3$y)
fd3.vp <- voronoi.polygons(fd3.vm)
fd3.vp
}
\keyword{ spatial }
