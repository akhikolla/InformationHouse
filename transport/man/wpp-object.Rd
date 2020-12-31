\name{wpp-object}
\alias{wpp-object}
\title{
  Class of Weighted Point Patterns
}
% maybe generalize later by a class wpp of weighted point patterns
\description{
  The class \code{"wpp"} represents discrete measures with positive mass at any of finitely many locations.
}
%\usage{
%FUNCTION_NAME(x, ...)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{x}{
%%     ~~Describe \code{x} here~~
%}
%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%}
%}
\details{
  Objects of class \code{"wpp"} may be created by the function 
  \code{\link{wpp}}, and are most commonly used as input to the function
  \code{\link{transport}}. There are methods \code{\link[=plot.wpp]{plot}}, \code{\link[=print.wpp]{print}} and
  \code{\link[=summary.wpp]{summary}} for this class.

  An object of class \code{"wpp"} contains the following elements:
  \tabular{ll}{
    \code{dimension} \tab
    the dimension of the Euclidean space in which the patterns live. Must be \eqn{\geq 2}{>=2}.
    \cr
    \code{N} \tab
    the total number of point.
    \cr
    \code{coordinates} \tab
    the coordinates of the points. An \code{N} \eqn{\times}{x} \code{dimension} matrix, where each row
    corresponds to a point.
    \cr
    \code{mass} \tab
    the masses at these points. A vector of length \code{N} of positive numbers.
    \cr
    \code{totmass} \tab
    the total mass of the point pattern.
  }  
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de} \cr
  Timo Wilm \email{timo.wilm@stud.uni-goettingen.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  Constructor function \code{\link{wpp}}.
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (x, ...) 
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{weight}
\keyword{weighted}
\keyword{point pattern}
\keyword{object}% __ONLY ONE__ keyword per line