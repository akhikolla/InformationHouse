\name{nearest.neighbours}
\alias{nearest.neighbours}
\title{
  Nearest neighbour structure for a data set
}
\description{
  This function can be used to generate nearest neighbour information
  for a set of 2D data points.
}
\usage{
nearest.neighbours(x, y)
}

\arguments{
  \item{x}{
    vector containing \eqn{x} ccordinates of points.
  }
  \item{y}{
    vector containing \eqn{x} ccordinates of points.
  }
}
\details{
  The C++ implementation of this function is used inside the
  \code{\link{locpoly}} and \code{\link{interp}} functions.
}
\value{
  A list with two components
  \item{index}{
    A matrix with one row per data point. Each row contains the indices
    of the nearest neigbours to the point associated with this row, 
    currently the point itself is also listed in the first row, so
    this matrix is of dimension \eqn{n} times \eqn{n} (will change to
    \eqn{n} times \eqn{n-1} later).
  }
  \item{dist}{
    A matrix containing the distances according to the neigbours listed
    in component \code{index}.
  }
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}


\seealso{
  \code{\link{convex.hull}}
}
\examples{
data(franke)
## use only a small subset
fd <- franke$ds1[1:5,]
nearest.neighbours(fd$x,fd$y)
}
\keyword{ utilities }

