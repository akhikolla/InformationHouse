\name{ibs}
\alias{ibs}
\title{
Integratal of a B-spline function
}
\description{
Calculate the integral of a B-spline function.
}
\usage{
ibs(x, knots, ord=4, coef = rep(1, length(knots) - ord))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Numerical value or vector. The value(s) at which to evaluate the
    integral of the B-spline; must be in the interval bewteen the
    smallest knot to the \code{(nknots-ord)}th largest knot.  
}
  \item{knots}{
    Numerical vector. The knot positions/sites of the B-spline function to be integrated.
}
  \item{ord}{
    An integer >=1. The order of the B-spline integrand function to be integrated. Equals degree plus 1.
}
  \item{coef}{
   A numerical vector. The coefficients (de Boor points) defining the B-spline integrand function.
}
}
\details{
  The function returns the integral(s) of the B-spline function
  specified  by knots \code{knots}, order \code{ord}, and coefficients
  \code{coef}, from the minimum knot position to each \code{x}
  value. The evaluation is based on a closed form expression of the
  integral in terms of higher order B-splines, given on page 128 of de
  Boor (2001).  
}
\value{
  A numerical equal to the integral(s).
}
\references{
  de Boor, C (2001) A Practical Guide to Splines. Revised Edition. Springer: New York.
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}

\seealso{
  \code{\link{bspline}}
}
\examples{
kns <- c(rep(0,4),1:4*0.2,rep(1,4))
co <- rnorm(length(kns)-3)
integrate(bspline,knots=kns,ord=3,coef=co,0,0.95)
integrate(function(x)bsbases(x,kns,3) \%*\% co,0,0.95)
ibs(0.95,kns,3,co)
}
\keyword{math}
