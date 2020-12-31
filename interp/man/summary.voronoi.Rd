\name{summary.voronoi}
\title{Return a summary of a voronoi object}
\usage{\method{summary}{voronoi}(object,...)
}
\alias{summary.voronoi}
\arguments{
  \item{object}{object of class \code{"voronoi"}}
  \item{...}{additional parameters for \code{summary}}
}
\description{
  Returns some information about \code{object}
}
\value{Object of class \code{"summary.voronoi"}.

  It contains the number of nodes (\code{nn}) and dummy nodes (\code{nd}).
}
\note{
  This function is meant as replacement for the function of same name in
  package \code{tripack} and should be fully backward compatible.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{voronoi}},\code{\link{voronoi.mosaic}},
  \code{\link{print.voronoi}}, \code{\link{plot.voronoi}},
  \code{\link{print.summary.voronoi}}.
}
\keyword{spatial}
