\name{plot.voronoi}
\title{Plot a voronoi object}
\usage{\method{plot}{voronoi}(x,add=FALSE,
                           xlim=c(min(x$tri$x)-
                             0.1*diff(range(x$tri$x)),
                             max(x$tri$x)+
                             0.1*diff(range(x$tri$x))),
                           ylim=c(min(x$tri$y)-
                             0.1*diff(range(x$tri$y)),
                             max(x$tri$y)+
                             0.1*diff(range(x$tri$y))),
                           all=FALSE,
                           do.points=TRUE,
                           main="Voronoi mosaic",
                           sub=deparse(substitute(x)),
                           isometric=TRUE,
                           ...)
}
\alias{plot.voronoi}
\arguments{
  \item{x}{object of class \code{"voronoi"}}
  \item{add}{logical, if \code{TRUE}, add to a current plot.}
  \item{xlim}{x plot ranges, by default modified to hide dummy 
        points outside of the plot}
  \item{ylim}{y plot ranges, by default modified to hide dummy 
        points outside of the plot}
  \item{all}{show all (including dummy points in the plot}
  \item{do.points}{logical, indicates if points should be plotted.}
  \item{main}{plot title}
  \item{sub}{plot subtitle}
  \item{isometric}{generate an isometric plot (default \code{TRUE})}
  \item{...}{additional plot parameters}
}
\description{Plots the mosaic \code{"x"}.
  Dashed lines are used for outer tiles of the mosaic.
}
\value{None
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{voronoi}}, \code{\link{print.voronoi}}, 
  \code{\link{summary.voronoi}}, \code{\link{plot.voronoi.polygons}} 
}
\examples{
data(franke)
tr <- tri.mesh(franke$ds3)
vr <- voronoi.mosaic(tr)
plot(tr)
plot(vr,add=TRUE)
}
\keyword{spatial}
