\name{plot.triSht}
\title{Plot a triangulation object}
\usage{\method{plot}{triSht}(x, add = FALSE, xlim = range(x$x),
  ylim = range(x$y), do.points = TRUE, do.labels = FALSE, isometric = TRUE,
  do.circumcircles = FALSE, segment.lty = "dashed", circle.lty =
  "dotted", ...)
}
\alias{plot.triSht}
\arguments{
  \item{x}{object of class \code{"triSht"}}
  \item{add}{logical, if \code{TRUE}, add to a current plot.}
  \item{do.points}{logical, indicates if points should be
    plotted. (default \code{TRUE})}
  \item{do.labels}{logical, indicates if points should be labelled.
    (default \code{FALSE})}
  \item{xlim,ylim}{x/y ranges for plot}
  \item{isometric}{generate an isometric plot (default \code{TRUE})}
  \item{do.circumcircles}{logical, indicates if circumcircles should be
    plotted (default \code{FALSE})}
  \item{segment.lty}{line type for triangulation segments}
  \item{circle.lty}{line type for circumcircles}
  \item{...}{additional plot parameters}
}
\description{plots the triangulation object \code{"x"}
}
\value{None
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}},
  \code{\link{summary.triSht}}
}
\examples{
## random points
plot(tri.mesh(rpois(100,lambda=20),rpois(100,lambda=20),duplicate="remove"))
## use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-10.78 & quakes[,1]>=-19.4 &
                     quakes[,2]<=182.29 & quakes[,2]>=165.77),]
quakes.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
plot(quakes.tri)
## use the whole quakes data set
## (will not work with standard memory settings, hence commented out)
## plot(tri.mesh(quakes$lon, quakes$lat, duplicate="remove"), do.points=F)
}
\keyword{spatial}
