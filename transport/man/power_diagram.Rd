\name{power_diagram}
\alias{power_diagram}
\alias{plot.power_diagram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Compute the Power Diagram of Weighted Sites in 2-Dimensional Space
}
\description{
   Compute the power diagram of weighted sites in 2-dimensional space.
}
\usage{
power_diagram(xi, eta, w, rect = NA)
\method{plot}{power_diagram}(x, weights=FALSE, add=FALSE, col=4, lwd=1.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi,eta,w}{
    vectors of equal length, where \code{xi}, \code{eta} are the coordinates of the sites and \code{w} are the corresponding weights.
}
  \item{rect}{
    vetor of length \code{4}. To get a finite representation of the power diagram, it will be intersected with the rectangle 
    \eqn{[rect[1],rect[2]] \times [rect[3],rect[4]]}{[rect[1],rect[3]]x[rect[2],rect[4]]}. Defaults to \code{c(min(xi),max(xi),min(eta),max(eta))}.
}
  \item{x}{
    a power diagram as returned from \code{\link{power_diagram}}.
}
  \item{weights}{
    logical. If \code{TRUE}, weights of non-redundant sites with non-negative weight are represented as circles whose radii
    are equal to the square roots of the corresponding weights.
}
  \item{add}{
    logical. Should the power diagram be plotted on top of current graphics?
}
  \item{col,lwd,\ldots}{
    further arguments graphic parameters used by \code{\link{plot.default}}.
}
}
\details{
  The function \code{power_diagram} implements an algorithm by Edelsbrunner and Shah (1996) which computes 
  regular triangulations and thus its dual representation, the power diagram. For point location, an algorithm 
  devised by Devillers (2002) is used.
}
\references{
  H. Edelsbrunner, N. R. Shah (1996), Incremental Topological Flipping Works for Regular Triangulations, Algorithmica 15, 223--241.

  O. Devillers (2002), The Delaunay Hierarchy, International Journal of Foundations of Computer Science 13, 163--180.
}
\author{
  Björn Bähre \email{bjobae@gmail.com} \cr
  (slightly modified by Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de})
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
    xi <- runif(100)
    eta <- runif(100)
    w <- runif(100,0,0.005)
    x <- power_diagram(xi,eta,w,rect=c(0,1,0,1))
    plot(x,weights=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
