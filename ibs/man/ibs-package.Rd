\name{ibs-package}
\alias{ibs-package}
\docType{package}
\title{
  A short title line describing what the package does
}
\description{
  Evaluate the integral of a B-spline function, or the B-spline function
  itself, or the B-spline basis functions.  
}
\details{
  The function \code{bsbases} calculates the values of the B-spline basis
  functions defined by knots \code{knots} and order \code{ord} at sites
  \code{x}. The function \code{bspline} evaluates the B-spline function
  defined by knots \code{knots}, order \code{ord}, and coefficients
  \code{coef} at sites \code{x}. The function \code{ibs} calculates the
  integrals of the B-spline function defined by knots \code{knots},
  order \code{ord}, and coefficients \code{coef} from the smallest knot
  to each of the values in \code{x}.
}
\author{
Creator: Feng Chen

Maintainer: Feng Chen <feng.chen@unsw.edu.au>
}
\references{
    de Boor, C. (2001) A Practical Guide to Splines. Revised Edition. Springer: New York.
}
\keyword{ package }
\keyword{ Splines}
\seealso{
  \code{\link[splines]{splineDesign}}
}
