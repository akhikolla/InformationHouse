\name{wpp}
\alias{wpp}
\title{
  Constructor for the wpp Class
}
\description{
  Construct an object of class \code{"wpp"} from a matrix of points and a vector of masses. 
}
\usage{
  wpp(coordinates, mass)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{
    a matrix specifying the coordinates of the points. Each row corresponds to a point.
  }
  \item{mass}{
  	a vector of non-negative values specifying the masses at these points.
  }
}
\details{
  For more detailed explanations of the arguments and other components of the returned object of class \code{"wpp"}, see 
  \code{\link{wpp-object}}. 
	
  It is legitimate to assign mass 0 to individual points in the arguments. However, when constructing the \code{wpp}-object such points are deleted. The coordinates of the deleted points can still be accessed via the attribute \code{zeropoints}. 
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
  Description of \link[=pp-object]{pp objects}.
}
\examples{
m <- matrix(c(1,1,2,2,3,1,4,2),4,2)
a <- pp(m)
print(a)
print.default(a)

\dontrun{  
  plot(a)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{weighted point pattern}
\keyword{object}% __ONLY ONE__ keyword per line