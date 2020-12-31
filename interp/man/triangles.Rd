\name{triangles}
\alias{triangles}
\title{Extract a list of triangles from a triangulation object}
\usage{triangles(tri.obj)
}
\arguments{
  \item{tri.obj}{object of class \code{\link{triSht}}}
}
\description{
  This function extracts a list of triangles
  from an triangulation object created by \code{tri.mesh}.
}
\details{
  The vertices in the returned matrix (let's denote it with
  \code{retval}) are ordered counterclockwise.  The columns \code{tr}\eqn{x}
  and \code{arc}\eqn{x}, \eqn{x=1,2,3} index the triangle and arc, respectively,
  which are opposite (not shared by) node \code{node}\eqn{x}, with
  \code{tri}\eqn{x=0} if \code{arc}\eqn{x} indexes a boundary arc.  Vertex indexes
  range from 1 to \eqn{n}, the number of nodes, triangle indexes from 0
  to \eqn{nt}, and arc indexes from 1 to \eqn{na = nt+n-1}.
}

\value{
  A matrix with columns \code{node1}, \code{node2}, \code{node3},
  representing the vertex nodal indexes,
  \code{tr1}, \code{tr2}, \code{tr3}, representing neighboring triangle
  indexes and \code{arc1}, \code{arc2}, \code{arc3} reresenting arc
  indexes.

  Each row represents one triangle.
}
%\references{
%}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}},
  \code{\link{plot.triSht}}, \code{\link{summary.triSht}},
  \code{\link{triangles}}
}
\examples{
# use the smallest Franke data set
data(franke)
fr3.tr<-tri.mesh(franke$ds3$x, franke$ds3$y)
triangles(fr3.tr)
}
\keyword{spatial}
