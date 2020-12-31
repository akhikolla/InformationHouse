\name{interp-package}
\alias{interp-package}
\docType{package}
\title{
  Interpolation of data
}
\description{
  Interpolation of \eqn{z} values given regular or irregular gridded
  data sets containing coordinates \eqn{(x_i,y_i)} and function values
  \eqn{z_i} is (will be) available through this package. As this
  interpolation is (for the irregular gridded data case) based on
  trianglation of the data locations also triangulation functions are
  implemented. Moreover the (not yet finished) spline interpolation
  needs estimators for partial derivates, these are also made available
  to the end user for direct use. 
}
\details{
  The interpolation use can be divided by the used method into piecewise
  linear (finished in 1_0.27) and spline  (not yet finished) interpolation and by
  input and output settings into  gridded and pointwise setups.
}
\note{
  This package is a FOSS replacement for the ACM licensed packages
  \code{akima} and \code{tripack}. The function calls are backward
  compatible.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>


  Maintainer: Albrecht Gebhardt <albrecht.gebhardt@aau.at>
}
%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
\keyword{ package }
\seealso{
  \code{\link{interp}}, \code{\link{tri.mesh}},
  \code{\link{voronoi.mosaic}}, \code{\link{locpoly}}
}
%\examples{
%  \dontrun{
%     ## Optional simple examples of the most important functions
%     ## These can be in \dontrun{} and \donttest{} blocks.   
%  }
%}
