\name{aha}
\alias{aha}
\alias{transport_apply}
\alias{transport_error}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Solve Transportation Problem by Aurenhammer--Hoffmann--Aronov Method
}
\description{
   Solve transportation problem by Aurenhammer--Hoffmann--Aronov Method.
}
\usage{
aha(a, b, nscales = 1, scmult = 2, factr = 1e+05, maxit = 10000, powerdiag=FALSE,
    wasser = FALSE, wasser.spt = NA, approx=FALSE, ...)
transport_apply(a, tplan)
transport_error(a, b, tplan)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
    an \eqn{m \times n}{m x n} matrix. \code{a} is treated as a measure on \eqn{[0,m] \times [0,n]}{[0,m] x [0,n]} with constant density on each unit square \eqn{[i,i+1) \times [j,j+1)}{[i,i+1)x[j,j+1)}.
}
  \item{b}{
    either a matrix such that \code{dim(a)==dim(b)} and \code{sum(a)==sum(b)} or a list of three vectors of equal length, named \code{x}, \code{y} and \code{mass} such that \code{sum(a)==sum(b$mass)},
    representing a discrete measure on \eqn{[0,m]) \times [0,n]}{[0,m]x[0,n]}.
}
  \item{tplan}{
    a transference plan from a (to b), typically an optimal transference plan obtained by a call to \code{aha}.  	
}
  \item{nscales, scmult}{
    the number of scales to use for the multiscale approach (the default is \eqn{1} meaning no multiscale approach), and the factor
    by which the number of pixels in each dimension is multiplied to get from a coarser to the next finer scale.
}
  \item{factr, maxit}{
    parameters passed to the underlying L-BFGS-B algorithm (via the argument \code{control} in the R-function \code{\link{optim}}).
}
  \item{powerdiag}{
    logical. Instead of an optimal transference plan, should the parameters for the optimal power diagram be returned? 
}
  \item{wasser}{
    logical. Instead of an optimal transference plan, should only the \eqn{L_2}{L2}-Wasserstein-distance between \code{a} and \code{b} be returned?
}
  \item{wasser.spt}{
  	the number of support points used to approximate the discrete measure \code{b}. Defaults to \code{NA} meaning the full set of support
  	points of \code{b} is used. If this argument is not \code{NA}, \code{wasser} is set to \code{TRUE}.
  }
  \item{approx}{
  	logical. If \code{TRUE}, an approximation to the objective function is used during optimization.
  }
  \item{\dots}{
    further arguments passed to \code{\link{optim}} via its argument \code{control}.
}
}
\details{
  The function \code{aha} implements the algorithm by Aurenhammer, Hoffmann and Aronov (1998) for finding optimal transference plans in terms
  of the squared Euclidean distance in two dimensions. It follows the more detailed description given in Mérigot (2011) and also implements
  the multiscale version presented in the latter paper.
  
  The functions \code{transport_apply} and \code{transport_error} serve for checking the accuracy of the transference plan obtained by \code{aha}.
  Since this transference plan is obtained by continuous optimization it will not transport exactly to the measure \code{b}, but to the measure
  \code{transport_apply(a, tplan)}. By \code{transport_error(a, b, tplan)} the sum of absolut errors between the transported \code{a}-measure and the \code{b}-measure is obtained.
}
\value{
If \code{powerdiag} and \code{wasser} are both \code{FALSE}, a data frame with columns \code{from}, \code{to} and \code{mass}, which specify from which knot to which other knot what amount of mass is sent in the optimal transference plan. Knots are given as indices in terms of the usual column major enumeration of the matrices \code{a} and \code{b}. There are \code{\link[=plot.pgrid]{plot}} methods for the classes \code{\link{pgrid}} and \code{\link{pp}}, which can plot this solution.

If \code{powerdiag} is TRUE and \code{wasser} is \code{FALSE}, a list with components \code{xi}, \code{eta}, \code{w} and \code{rect}, which specify the parameters for the optimal power diagram in the same format as needed for the function \code{\link{power_diagram}}. Note that rect is always \code{c(0,m,0,n)}. Since version 0.10-0 the list has a further component \code{wasser.dist} containing the Wasserstein distance.

If \code{wasser} is \code{TRUE}, a data frame with columns \code{wasser.dist} and \code{error.bound} of length one, where \code{error.bound} gives a bound on the absolute error in the Wasserstein distance due to approximating the measure \code{b} by a measure on a smaller number of support points.
}

\references{
  F. Aurenhammer, F. Hoffmann and B. Aronov (1998). Minkowski-type theorems and least-squares clustering. Algorithmica 20(1), 61--76.

  Q. Mérigot (2011). A multiscale approach to optimal transport. Eurographics Symposium on Geometry Processing 30(5), 1583--1592.
}
\author{
  Björn Bähre \email{bjobae@gmail.com} \cr
  (slightly modified by Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de})
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{transport}}, which is a convenient wrapper function for various optimal transportation algorithms.
}
\examples{
res <- aha(random32a$mass, random32b$mass)
plot(random32a, random32b, res, lwd=0.75)

aha(random64a$mass, random64b$mass, nscales=3, scmult=5, wasser.spt=512, approx=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
