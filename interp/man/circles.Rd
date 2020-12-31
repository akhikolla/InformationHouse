\name{circles}
\alias{circles}
\title{Add circles to a plot}
\description{
 This function plots circles at given locations with given radii.
}
\usage{
circles(x, y, r, ...)
}
\arguments{
  \item{x}{ vector of x coordinates }
  \item{y}{ vector of y coordinates }
  \item{r}{ vactor of radii }
  \item{\dots}{ additional graphic parameters will be passed through }
}
\note{
  This function needs a previous plot where it adds the circles.

  This function was earlier used in package \code{tripack}.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{
  \code{\link{lines}},
  \code{\link{points}}
}
\examples{
x<-rnorm(10)
y<-rnorm(10)
r<-runif(10,0,0.5)
plot(x,y, xlim=c(-3,3), ylim=c(-3,3), pch="+")
circles(x,y,r)
}
\keyword{ aplot }

