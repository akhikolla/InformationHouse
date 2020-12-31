\name{bsbases}
\alias{bsbases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  B-spline bases
}
\description{
  Evaluate the B-spline basis function at the \code{x} value(s). 
}
\usage{
bsbases(x, knots, ord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Numerical value or vector. The value(s) at which to evaluate the B
    spline bases.
}
  \item{knots}{
    Numerical vector. The knot positions/sites of the B-spline bases.
}
  \item{ord}{
    An integer >=1. The order of the B-spline bases. Equals degree plus 1.
}
}

\value{
  A numerical matrix of \code{length(x)} rows and
  \code{length(knots)-ord} columns.
}
\references{
  de Boor, C. (2001) A Practical Guide to Splines. Revised Edition. Springer: New York.
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\note{
  In contrast to the implementation based on the \code{splineDesign}
  function from the \code{splines} package, this version of the B-spline
  basis functions are left-continuous at the rightmost knot.   
}

\seealso{
    \code{\link[splines]{splineDesign}}
}
\examples{
kns <- c(rep(0,4),1:4*0.2,rep(1,4))
round(bsbases(0:5/5,kns,2),10)
round(splines::splineDesign(kns,0:5/5,2),10)
}
\keyword{ math }

