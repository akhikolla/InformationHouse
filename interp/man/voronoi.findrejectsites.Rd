\name{voronoi.findrejectsites}
\title{Find the Voronoi sites at the border of the region (to be rejected).}
\author{S. J. Eglen}
\usage{voronoi.findrejectsites(voronoi.obj, xmin, xmax, ymin, ymax)
}
\alias{voronoi.findrejectsites}
\arguments{
  \item{voronoi.obj}{object of class \code{"voronoi"}}
  \item{xmin}{minimum x-coordinate of sites in the region}
  \item{xmax}{maximum x-coordinate of sites in the region}
  \item{ymin}{minimum y-coordinate of sites in the region}
  \item{ymax}{maximum y-coordinate of sites in the region}

}
\description{Find the sites in the Voronoi tesselation that lie at the
  edge of the region.  A site is at the edge if any of the vertices of
  its Voronoi polygon lie outside the rectangle with corners (xmin,ymin)
  and (xmax,ymax).
}
\value{A logical vector of the same length as the number of sites.  If
  the site is a reject, the corresponding element of the vector is set
  to TRUE.}
\seealso{
  \code{\link{voronoi.polygons}}
}
\keyword{spatial}

