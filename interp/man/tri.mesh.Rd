\name{tri.mesh}
\alias{tri.mesh}
\title{
  Delaunay triangulation
}
\description{
  This function generates a Delaunay triangulation of arbitrarily distributed
  points in the plane. The resulting object can be printed or plotted,
  some additional functions can extract details from it like the list of
  triangles, arcs or the convex hull.
}
\usage{
tri.mesh(x, y = NULL, duplicate = "error")
}
\arguments{
  \item{x}{
    vector containing \eqn{x} coordinates of the data. If \code{y} is missing
    \code{x} should be a list or dataframe with two components \code{x}
    and \code{y}.
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
  This function creates a Delaunay triangulation of a set of
  arbitrarily distributed points in the plane referred to as
  nodes.

  The Delaunay triangulation is defined as a set of triangles with the
  following five properties:
  \enumerate{
    \item The triangle vertices are nodes.
    \item No triangle contains a node other than its vertices.
    \item The interiors of the triangles are pairwise disjoint.
    \item The union of triangles is the convex hull of the set
    of nodes (the smallest convex set which contains
    the nodes).
    \item The interior of the circumcircle of each triangle
    contains no node.
    }
    The first four properties define a triangulation, and the
    last property results in a triangulation which is as close
    as possible to equiangular in a certain sense and which is
    uniquely defined unless four or more nodes lie on a common
    circle.  This property makes the triangulation well-suited
    for solving closest point problems and for triangle-based
    interpolation.

    This triangulation is based on the s-hull algorithm by David
    Sinclair. It consist of two steps:
    \enumerate{
      \item{
        Create an initial non-overlapping triangulation from the
        radially sorted nodes (w.r.t to an arbitrary first
        node). Starting from a first triangle built from the first node
        and its nearest neigbours this is done by adding triangles from the
        next node (in the sense of distance to the first node) to the hull
        of the actual triangulation visible from this node (sweep hull
        step).
      }
      \item{
        Apply triange flipping to each pair of triangles sharing a
        border until condition 5 holds (Cline-Renka test).
      }
    }
    This algorithm has complexicity \eqn{O(n*log(n))}.
}
\value{
  an object of class \code{"triSht"}, see \code{\link{triSht}}.
}
\references{
  B. Delaunay, Sur la sphere vide. A la memoire de Georges Voronoi,
  Bulletin de l'Academie des Sciences de l'URSS. Classe des sciences
  mathematiques et na, 1934, no. 6, p. 793--800

  D. A. Sinclair, S-Hull: A Fast Radial Sweep-Hull Routine for Delaunay
  Triangulation. https://arxiv.org/pdf/1604.01428.pdf, 2016.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\note{
  This function is meant as a replacement for
  \code{\link[tripack]{tri.mesh}} from package \code{tripack}.
  Please note that the underlying algorithm changed from Renka's method
  to Sinclair's sweep hull method. Delaunay triangulations are unique if
  no four or more points exist which share the same
  circumcircle. Otherwise several solutions are available and different
  algorithms will give different results. This especially holds for
  regular grids, where in the case of rectangular gridded points each
  grid cell can be triangulated in two different ways.

  The arguments are backward compatible, but the returned object
  is not compatible with package \code{tripack} (it
  provides a \code{\link[tripack]{tri}} object type)! But you
  can apply methods with same names to the object returned in package
  \code{interp} which is of type \code{\link{triSht}}, so you can reuse
  your old code but you cannot reuse your old saved workspace.
}


\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}},
  \code{\link{summary.triSht}}, \code{\link{triangles}},
  \code{\link{convex.hull}}, \code{\link{arcs}}.
}
\examples{
## use Frankes datasets:
data(franke)
tr1 <- tri.mesh(franke$ds3$x, franke$ds3$y)
tr1
tr2 <- tri.mesh(franke$ds2)
summary(tr2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

