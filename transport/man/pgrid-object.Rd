\name{pgrid-object}
\alias{pgrid-object}
\title{
  Class of Pixel Grids
}
\description{
  The class \code{"pgrid"} (for pixel grid) represents regular quantizations of measures on
  (bounded subsets of) \eqn{R^d}. Currently only square quantizations of measures on a rectangles
  are supported, which in 2-d can be thought of as grey scale images. 
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
  Objects of class \code{"pgrid"} can be created by the function 
  \code{\link{pgrid}}, and are most commonly used as input to the function
  \code{\link{transport}}. There are methods \code{\link[=plot.pgrid]{plot}}, \code{\link[=print.pgrid]{print}} and
  \code{\link[=summary.pgrid]{summary}} for this class.

  An object of class \code{"pgrid"} contains the following elements:
  \tabular{ll}{
    \code{structure} \tab
    the structure of the grid. Currently only \code{"square"} and \code{"rectangular"} are supported.
    \cr
    \code{dimension} \tab
    the dimension \eqn{d} of the space in which the grid is embedded. Must be \eqn{\geq 2}{>=2}.
    \cr
        \code{n} \tab
    the number of pixels along the various coordinates, a vector of length \code{dimension}.
    \cr
    \code{N} \tab
    the total number of pixels.
    \cr
    \code{boundary} \tab
    the outer boundary of the "picture" (i.e. of the support of the measure). A vector of \cr \tab
    length \code{2*dimension}, where the odd entries contain the left and the even entries \cr \tab
    contain the right endpoints of the various coordinates.
    \cr
    \code{gridtriple} \tab
    the rule for generating the pixel centres along the various coordinates. \cr \tab
    A \code{dim} by \code{3} matrix where each row is of the form \code{c(start, end, step)}.
    \cr
    \code{generators} \tab
    the pixel centres along the various coordinates. A list of length \code{dim} where the \cr \tab
    \code{i}-th element is a vector of length \code{n[i]}. 
    \cr
    \code{mass} \tab
    the array of masses in each pixel / at each pixel centre. In 2-d orientation corresponds \cr \tab
    to the standard orientation of images, see e.g. \code{\link{image}}. This means that \cr \tab
    pixels are arranged on coordinate axes in the order of their indices.
    \cr
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
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  Constructor function \code{\link{pgrid}}.
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
\keyword{grid}
\keyword{object}% __ONLY ONE__ keyword per line