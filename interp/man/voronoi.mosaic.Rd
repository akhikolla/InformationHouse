\name{voronoi.mosaic}
\alias{voronoi.mosaic}
\title{
  Voronoi mosaic
}
\description{
    This function creates a Voronoi mosaic out of a given set of
    arbitraryly located points in the plane. Each cell of a voronoi
    mosaic is associated with a data point and contains all points
    \eqn{(x,y)} closest to this data point.
}
\usage{
voronoi.mosaic(x, y = NULL, duplicate = "error")
}
\arguments{
  \item{x}{
    vector containing \eqn{x} coordinates of the data. If \code{y} is missing
    \code{x} should be a list or dataframe with two components \code{x}
    and \code{y}.

    \code{x} can also be an object of class \code{\link{triSht}} generated
    by \code{\link{tri.mesh}}. In this case the internal triangulation
    step can be skipped.
  }
  \item{y}{
    vector containing \eqn{y} coordinates of the data. Can be omitted if
    \code{x} is a list with two components \code{x} and \code{y}.
  }
  \item{duplicate}{
    flag indicating how to handle duplicate elements.
    Possible values are:
    \itemize{
      \item{
        \code{"error"} -- default,
      }
      \item{
        \code{"strip"} -- remove all duplicate points,
      }
      \item{
        \code{"remove"} -- leave one point of the duplicate points.
      }
    }
  }
}
\details{
  The function creates first a Delaunay triangulation (if not already
  given), extracts the circumcircle centers of these triangles, and then
  connects these points according to the neighbourhood relations between
  the triangles.
}
\value{
  An object of class \code{\link{voronoi}}.
}
\references{
G. Voronoi, Nouvelles applications des parametres continus a la theorie
    des formes quadratiques. Deuxieme memoire. Recherches sur les
    parallelloedres primitifs, Journal fuer die reine und angewandte
    Mathematik, 1908, vol 134, p. 198-287
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\note{
  This function is meant as a replacement for
  \code{\link[tripack]{voronoi.mosaic}} from package \code{tripack}.
  Please note that the underlying triangulation uses a
  different algorithm, see \code{\link{tri.mesh}}. Contrary to
  \code{tri.mesh} this should not affect the result for non
  unique triangulations e.g. on regular grids as the voronoi mosaic in
  this case will still be unique.

  The arguments are backward compatible, even the returned object should be
  compatible with functions from package \code{tripack}.
}

\seealso{
  \code{\link{voronoi}},\code{\link{voronoi.mosaic}},
  \code{\link{print.voronoi}}, \code{\link{plot.voronoi}}
}
\examples{
data(franke)
fd <- franke$ds3
vr <- voronoi.mosaic(fd$x, fd$y)
summary(vr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }% use one of  RShowDoc("KEYWORDS")

