\name{voronoi.area}
\title{Calculate area of Voronoi polygons}
\author{S. J. Eglen}
\usage{voronoi.area(voronoi.obj)
}
\alias{voronoi.area}
\arguments{
  \item{voronoi.obj}{object of class \code{"voronoi"}}
}
\description{Computes the area of each Voronoi polygon.
For some sites at the edge of the region, the Voronoi polygon is not
bounded, and so the area of those sites cannot be calculated, and hence
will be \code{NA}.
}
\value{A vector of polygon areas.}
\seealso{
  \code{\link{voronoi.mosaic}},\code{\link{voronoi.polygons}},
}
\keyword{spatial}

\examples{
data(franke)
fd3 <- franke$ds3
fd3.vm <- voronoi.mosaic(fd3$x,fd3$y)
fd3.vm.areas <- voronoi.area(fd3.vm)
plot(fd3.vm)
text(fd3$x, fd3$y, round(fd3.vm.areas,5))
}

