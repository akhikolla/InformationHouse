\name{tri.find}
\title{Locate a point in a triangulation}
\usage{
tri.find(tri.obj,x,y)
}
\alias{tri.find}
\arguments{
  \item{tri.obj}{an triangulation object of class \code{triSht}}
    \item{x}{x-coordinate of the point}
    \item{y}{y-coordinate of the point}
}
\description{
  This subroutine locates a point \eqn{P=(x,y)} relative to a
  triangulation  created by \code{tri.mesh}.  If \eqn{P} is
  contained in a triangle, the three vertex indexes are
  returned.  Otherwise, the indexes of the rightmost and
  leftmost visible boundary nodes are returned.
}
\value{
  A list with elements \code{i1},\code{i2},\code{i3} containing nodal
  indexes, in counterclockwise order, of the vertices of a triangle
  containing \eqn{P=(x,y)}. \code{bc} contains the barycentric coordinates 
  of \eqn{P} w.r.t. the found triangle.

  If \eqn{P} is not contained in the
  convex hull of the nodes this indices are 0 (\code{bc} is meaningless then).
  % CHECKME (maybe differs from tripack!):
  %, \code{i1} indexes the rightmost visible
  %boundary node, \code{i2} indexes the leftmost visible boundary node,
  %and \code{i3} = 0.  Rightmost and leftmost are defined from the
  %perspective of \eqn{P}, and a pair of points are visible from each other if
  %and only if the line segment joining them intersects no triangulation
  %arc.  If \eqn{P} and all of the nodes lie on a common line, then
  %\code{i1}=\code{i2}=\code{i3} = 0 on output.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}

\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}},
  \code{\link{summary.triSht}}, \code{\link{triangles}},
  \code{\link{convex.hull}}
}
\examples{
data(franke)
tr<-tri.mesh(franke$ds3$x,franke$ds3$y)
plot(tr)
pnt<-list(x=0.3,y=0.4)
triangle.with.pnt<-tri.find(tr,pnt$x,pnt$y)
attach(triangle.with.pnt)
lines(franke$ds3$x[c(i1,i2,i3,i1)],franke$ds3$y[c(i1,i2,i3,i1)],col="red")
points(pnt$x,pnt$y)
}
\keyword{spatial}
