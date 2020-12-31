\name{on.convex.hull}
\title{Determines if points are on or in the convex hull of a triangulation object}
\usage{
on.convex.hull(tri.obj, x, y, eps=1E-16)
in.convex.hull(tri.obj, x, y, eps=1E-16, strict=TRUE)
}
\alias{on.convex.hull}
\alias{in.convex.hull}
\arguments{
  \item{tri.obj}{object of class \code{\link{triSht}}}
  \item{x}{vector of \eqn{x}-coordinates of points to locate}
  \item{y}{vector of \eqn{y}-coordinates of points to locate}
  \item{eps}{accuracy for checking the condition}
  \item{strict}{logical, default \code{TRUE}. It indicates if the convex
    hull is treated as an open (\code{strict=TRUE}) or closed
    (\code{strict=FALSE}) set. (applies only to \code{in.convex.hull})}
}
\description{
Given a triangulation object \code{tri.obj} of \eqn{n} points in the plane, this
subroutine returns a logical vector indicating if the points
\eqn{(x_i,y_i)}
lay on or in the convex hull of \code{tri.obj}.

}
\value{
Logical vector.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}},
  \code{\link{summary.triSht}}, \code{\link{triangles}},
  \code{\link{convex.hull}}.
}
\examples{
# use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-10.78 & quakes[,1]>=-19.4 &
                     quakes[,2]<=182.29 & quakes[,2]>=165.77),]
q.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
on.convex.hull(q.tri,quakes.part$lon[1:20],quakes.part$lat[1:20])
# Check with part of data set:
# Note that points on the hull (see above) get marked FALSE below:
in.convex.hull(q.tri,quakes.part$lon[1:20],quakes.part$lat[1:20])
# If points both on the hull and in the interior of the hull are meant 
# disable strict mode:
in.convex.hull(q.tri,quakes.part$lon[1:20],quakes.part$lat[1:20],strict=FALSE)
# something completely outside:
in.convex.hull(q.tri,c(170,180),c(-20,-10))
}
\keyword{spatial}
