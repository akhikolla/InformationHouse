\name{wasserstein1d}
\alias{wasserstein1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the Wasserstein Distance Between Two Univariate Samples
}
\description{
  Given two vectors \code{a} and \code{b}, compute the Wasserstein distance of
  order \code{p} between their empirical distributions.
}
\usage{
wasserstein1d(a, b, p = 1, wa = NULL, wb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, b}{
    two vectors.
  }
  \item{p}{
    a positive number. The order of the Wasserstein distance.
  }
  \item{wa, wb}{
    optional vectors of weights for \code{a} and \code{b}. 
  }
}
\details{
  The Wasserstein distance of order \code{p} is defined as the \code{p}-th root of the total cost incurred when transporting a pile of mass into another pile of mass in an optimal way, where the cost of transporting a unit of mass from \eqn{x} to \eqn{y} is given as the \code{p}-th power \eqn{\|x-y\|^p}{||x-y||^p} of the Euclidean distance. 
  
  In the present function the vector \code{a} represents the locations on the real line of \eqn{m} deposits of mass \eqn{1/m} and the vector \code{b} the locations of \eqn{n} deposits of mass \eqn{1/n}. If the user specifies weights \code{wa} and \code{wb}, these default masses are replaced by \code{wa/sum(wa)} and \code{wb/sum(wb)}, respectively.
  
  In terms of the empirical distribution function \eqn{F(t) = \sum_{i=1}^m w^{(a)}_i 1\{a_i \leq t\}}{F^{(a)}(t) = sum_{i=1}^m w^{(a)}_i 1{a_i \leq t}} of locations \eqn{a_i} with normalized weights \eqn{w^{(a)}_i}, and the corresponding function \eqn{G(t) = \sum_{j=1}^n w^{(b)}_j 1\{b_j \leq t\}}{G(t) = sum_{j=1}^n w^{(b)}_j 1{b_j \leq t}} for \code{b}, the Wasserstein distance in 1-d is given as
  \deqn{W_p(F,G) = \left(\int_0^1 |F^{-1}(u)-G^{-1}(u)|^p \; du \right)^{1/p},}{W_p(F,G) = ( integral_0^1 |F^{-1}(u)-G^{-1}(u)|^p du )^{1/p},}
  where \eqn{F^{-1}} and \eqn{G^{-1}} are generalized inverses. If \eqn{p=1}, we also have \deqn{W_1(F,G) = \int_{-\infty}^{\infty} |F(x)-G(x)| \; dx.}{W_1(F,G) = ( integral_{-Inf}^{Inf} |F(x)-G(x)| dx).}
}
\value{
  A single number, the Wasserstein distance for the specified data.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{wasserstein}}
}
\examples{
  x <- rnorm(200)
  y <- rnorm(150,2)
  wasserstein1d(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mallows}
\keyword{earth mover}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
