\name{identify.triSht}
\title{Identify points in a triangulation plot}
\usage{\method{identify}{triSht}(x,...)
}
\alias{identify.triSht}
\arguments{
  \item{x}{object of class \code{\link{triSht}}}
  \item{...}{additional paramters for \code{identify}}
}
\description{Identify points in a plot of \code{"x"} with its
  coordinates. The plot of \code{"x"} must be generated with \code{plot.tri}.
}
\value{an integer vector containing the indexes of the identified points.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
    \code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}}, \code{\link{summary.triSht}}
}
\examples{
\dontrun{
data(franke)
tr <- tri.mesh(franke$ds3$x, franke$ds3$y)
plot(tr)
identify(tr)
}
}
\keyword{spatial}
