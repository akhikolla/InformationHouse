\name{wasserstein}
\alias{wasserstein}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the Wasserstein Distance Between Two Objects
}
\description{
  Given two objects \code{a} and \code{b} that specify measures in \eqn{R^d}, compute the Wasserstein distance of
  order \code{p} between the objects.
}
\usage{
wasserstein(a, b, p=1, tplan=NULL, costm=NULL, prob=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, b}{
    two objects that describe mass distributions in \eqn{R^d}. Either both of class \code{\link{pgrid}} or \code{\link{pp}}
    or \code{\link{wpp}} or \code{numeric}. For the first three the dimension \eqn{d} of the structures must be at least 2; see function \code{\link{wasserstein1d}} for \eqn{d=1}.
}
  \item{p}{
    the power \eqn{\geq 1}{>=1} to which the Euclidean distance between points is taken in order to
    compute transportation costs.
}
  \item{tplan}{
    an optional transference plan in the format returned by the function \code{\link{transport}}. If \code{NULL} an optimal
    transference plan based on \code{a}, \code{b} and \code{p} is computed by a call to \code{\link{transport}}.
}
  \item{costm}{
    the matrix of costs between the support points of the measures. Ignored unless \code{a} and \code{b} are numeric vectors.
}
  \item{prob}{
  	logical. Should the objects \code{a}, \code{b} be interpreted as probability measures, i.e. their total mass be normalized to 1?
  }
  \item{...}{
  	further parameters passed to \code{\link{transport}} if \code{tplan} is \code{NULL}.
  }
}
\details{
  The Wasserstein distance of order \code{p} is defined as the \code{p}-th root of the total cost incurred when transporting measure \code{a}
  to measure \code{b} in an optimal way, where the cost of transporting a unit of mass from \eqn{x} to \eqn{y} is given as the \code{p}-th 
  power \eqn{\|x-y\|^p}{||x-y||^p} of the Euclidean distance. 
	
  If \code{tplan} is supplied by the user, no checks are performed whether it is optimal for the given problem. So this
  function may be used to compare different (maybe suboptimal) transference plans with regard to their total costs.
  
  For further details on the algorithms used, see help of \code{\link{transport}}. 
}
\value{
  A single number, the Wasserstein distance for the specified problem.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=plot.pgrid]{plot}}, \code{\link{transport}}, \code{\link{wasserstein1d}}
}
\examples{
#
# example for class 'pgrid'
#
wasserstein(random32a,random32b,p=1)
res <- transport(random32a,random32b,p=2)
wasserstein(random32a,random32b,p=1,res)
# is larger than above:
# the optimal transport for p=2 is not optimal for p=1

#
# example for class 'pp'
#
set.seed(27)
x <- pp(matrix(runif(500),250,2))
y <- pp(matrix(runif(500),250,2))
wasserstein(x,y,p=1)
wasserstein(x,y,p=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mallows}
\keyword{earth mover}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
