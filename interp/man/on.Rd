\name{on}
\alias{on}

\alias{left}

\title{
  Determines if a point is on or left of the vector described by two other points.
}
\description{
  A simple test function to determine the position of one (or more)
  points relative to a vector spanned by two points.
}
\usage{
on(x1, y1, x2, y2, x0, y0, eps = 1e-16)
left(x1, y1, x2, y2, x0, y0, eps = 1e-16)
}
\arguments{
  \item{x1}{
    \code{x} coordinate of first point determinig the vector.
  }
  \item{y1}{
    \code{y} coordinate of first point determinig the vector.
  }
  \item{x2}{
    \code{x} coordinate of second point determinig the vector.
  }
  \item{y2}{
    \code{y} coordinate of second point determinig the vector.
  }
  \item{x0}{
    vector of \code{x} coordinates to locate relative to the vector
    \eqn{(x_2-x_1, y_2-y_1)}.
  }
  \item{y0}{
    vector of \code{x} coordinates to locate relative to the vector
    \eqn{(x_2-x_1, y_2-y_1)}.
  }
  \item{eps}{
    tolerance for checking if \eqn{x_0,y_0} is on or left of \eqn{(x_2-x_1,
      y_2-y_1)}, defaults to \eqn{10^{-16}}.
  }
}
\value{
  logical vector with the results of the test.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}

\seealso{
  \code{\link{in.convex.hull}},   \code{\link{on.convex.hull}}.
}
\examples{
y <- x <- c(0,1)
## should be TRUE
on(x[1],y[1],x[2],y[2],0.5,0.5)
## note the default setting of eps leading to
on(x[1],y[1],x[2],y[2],0.5,0.50000000000000001)
## also be TRUE

## should be TRUE
left(x[1],y[1],x[2],y[2],0.5,0.6)
## note the default setting of eps leading to
left(x[1],y[1],x[2],y[2],0.5,0.50000000000000001)
## already resulting to FALSE

}
\keyword{ utilities }

