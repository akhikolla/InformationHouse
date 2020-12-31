\name{arcs}
\alias{arcs}
\title{
  Extract a list of arcs from a triangulation object.
}
\description{
  This function extracts a list of arcs from a triangulation object
  created by \code{tri.mesh}.
}
\usage{
arcs(tri.obj)
}
\arguments{
  \item{tri.obj}{
    object of class \code{\link{triSht}}
  }
}
\details{
  This function acesses the \code{arcs} component of a triangulation
  object returned by \code{\link{tri.mesh}} and extracts the arcs
  contained in this triangulation. This is e.g. used for plotting.
}
\value{
  A matrix with two columns \code{"from"} and \code{"to"} containing the
  indices of points connected by the arc with the corresponding row index.

}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{triangles}}, \code{\link{area}}
}
\examples{
data(franke)
tr <- tri.mesh(franke$ds3)
arcs(tr)
}
\keyword{ spatial }
\keyword{ dplot }

