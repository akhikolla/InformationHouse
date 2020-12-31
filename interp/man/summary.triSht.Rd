\name{summary.triSht}
\title{Return a summary of a triangulation object}
\usage{\method{summary}{triSht}(object,...)
}
\alias{summary.triSht}
\arguments{
  \item{object}{object of class \code{"triSht"}}
  \item{...}{additional paramters for \code{summary}}
}
\description{
  Returns some information (number of nodes, triangles, arcs)
  about \code{object}.
}
\value{An object of class \code{"summary.triSht"}, to be printed by
  \code{\link{print.summary.triSht}}.

  It contains the number of nodes (\code{n}), of arcs (\code{na}), of
  boundary nodes (\code{nb}) and triangles (\code{nt}).
}
\note{
  This function is meant as replacement for the function of same name in
  package \code{tripack}.

  The only difference is that no constraints are possible with
  \code{triSht} objects of package \code{interp}.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{triSht}}, \code{\link{print.triSht}}, \code{\link{plot.triSht}},
  \code{\link{print.summary.triSht}}.
}
\keyword{spatial}
