\name{bspline}
\alias{bspline}
\title{
  B-spline function
}
\description{
  Evaluate the value of the B-spline function with knots \code{knots},
  order \code{ord}, and coefficients \code{coef} at value(s) \code{x}. 
}
\usage{
bspline(x, knots, ord, coef = rep(1, length(knots) - ord))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Numerical value or vector. The value(s) at which to evaluate the B-spline.
}
  \item{knots}{
    Numerical vector. The knot positions/sites of the B-spline function.
}
  \item{ord}{
    An integer >=1. The order of the B-spline function. Equals degree plus 1.
}
  \item{coef}{
   A numerical vector. The coefficients (de Boor points) defining the B-spline function.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A scalar or a vector of length equal to that of \code{x}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  function is left-continuous at the rightmost knot site.   
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[splines]{splineDesign}}
}
\examples{
kns <- c(rep(0,4),1:4*0.2,rep(1,4))
bspline(0:10/10,kns,3,rep(1,length(kns)-3))
splines::splineDesign(kns,0:10/10,3)\%*\%rep(1,length(kns)-3)
}
\keyword{math}