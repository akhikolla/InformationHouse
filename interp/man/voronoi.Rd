\name{voronoi}
\title{Voronoi object}
\alias{voronoi}
\arguments{
  \item{x,y}{x and y coordinates of nodes of the voronoi mosaic. Each node is a
    circumcircle center of some triangle from the Delaunay triangulation.}
  \item{node}{logical vector, indicating real nodes of the voronoi
    mosaic. These nodes are the centers of circumcircles of triangles with
    positive area of the delaunay triangulation.

    If \code{node[i]}=\code{FALSE}, (\code{c[i]},\code{x[i]}) belongs to a
    triangle with area 0.}
  \item{n1,n2,n3}{indices of neighbour nodes. Negative indices indicate
    dummy points as neighbours.}
  \item{tri}{triangulation object, see \code{\link{triSht}}.}
  \item{area}{area of triangle \eqn{i}.
  }
  \item{ratio}{aspect ratio (inscribed radius/circumradius) of triangle
    \eqn{i}.}
  \item{radius}{circumradius of triangle i.}
  \item{dummy.x,dummy.y}{x and y coordinates of dummy points. They are
    used for plotting of unbounded tiles.}
}
\description{
  A \code{voronoi} object is created with \code{\link{voronoi.mosaic}}
}
\note{
  This version of \code{voronoi} object is generated from the
  \code{\link{tri.mesh}} function from package \code{interp}. That's the only
  difference to \code{voronoi} objects generated with package
  \code{tripack}.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{voronoi.mosaic}},\code{\link{plot.voronoi}}
}
\keyword{spatial}

