\name{pgrid}
\alias{pgrid}
\title{
  Constructor for the pgrid Class
}
\description{
  Construct an object of class \code{"pgrid"} from a matrix or a higher-dimensional array. 
}
\usage{
  pgrid(mass, boundary, gridtriple, generator, structure) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{
    a matrix or higher-dimensional array specifing the masses in each pixel / at each pixel centre.
}
  \item{boundary, gridtriple, generator}{
    arguments specifying the positions of the pixels. At most one of these can be specified.
}
  \item{structure}{
    optional character string specifying the structure of the grid.
    Currently only \code{"square"} and \code{"rectangular"} make sense, and are derrived automatically
    from the dimensions of \code{mass}.
}
}
\details{
  For more detailed explanations of the arguments and other components of the derived object of class \code{"pgrid"}, see 
  \code{\link{pgrid-object}}. 
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
  Description of \link[=pgrid-object]{pgrid objects}.
}
\examples{
m <- matrix(1:20, 4, 5)
a <- pgrid(m)
print(a)
print.default(a)

\dontrun{  
  plot(a, rot=TRUE)}
}
\keyword{grid}
\keyword{object}% __ONLY ONE__ keyword per line