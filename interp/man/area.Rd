\name{area}
\alias{area}
\title{
  Extract a list of triangle areas from a triangulation object.
}
\description{
  This function returns a list containing the areas of each triangle of
  a  triangulation object created by \code{tri.mesh}.
}
\usage{
area(tri.obj)
}
\arguments{
  \item{tri.obj}{
    object of class \code{\link{triSht}}
  }
}
\details{
  This function acesses the \code{cclist} component of a triangulation
  object returned by \code{\link{tri.mesh}} and extracts the areas
  of the triangles contained in this triangulation.
}
\value{
  A vector containing the area values.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{triangles}}, \code{\link{arcs}}
}
\examples{
data(franke)
tr <- tri.mesh(franke$ds3)
area(tr)
}
\keyword{ spatial }
\keyword{ utilities }

