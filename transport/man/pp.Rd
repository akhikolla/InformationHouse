\name{pp}
\alias{pp}
\title{
  Constructor for the pp Class
}
\description{
  Construct an object of class \code{"pp"} from a matrix. 
}
\usage{
  pp(coordinates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{
    a matrix specifying the coordinates of the points. Each row corresponds to a point.
  }
}
\details{
  For more detailed explanations of the arguments and other components of the derived object of class \code{"pp"}, see 
  \code{\link{pp-object}}. 
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
\keyword{point pattern}
\keyword{object}% __ONLY ONE__ keyword per line