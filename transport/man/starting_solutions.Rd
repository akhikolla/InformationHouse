\name{starting solutions}
\alias{northwestcorner}
\alias{russell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute starting solution for the transportation problem
}
\description{
  Compute a feasible transference plan between two mass vectors.
}
\usage{
northwestcorner(a, b)
russell(a, b, costm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a,b}{
    Two numeric vectors (typically containing natural numbers) of length \eqn{m} and \eqn{n}, describing mass distributions.
}
  \item{costm}{
    A \eqn{m} by \eqn{n} matrix of costs for moving one unit of mass.
}
}
\value{
A list whose components are \eqn{m} by \eqn{n} matrices, viz.
  \item{assignment}{containing as \eqn{(i,j)}-th entry the mass assigned from origin \eqn{i} to destination \eqn{j};}
  \item{basis}{containing as \eqn{(i,j)}-th entry a \eqn{1} if it is a basic entry and a \eqn{0} otherwise.}
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}

\section{Warnings}{The current implementations are in R. Computations may be slow for larger vectors \code{a} and \code{b}.\cr
	The computed starting solution may be degenerate, i.e. there may be basic entries where zero mass is assigned.}

\seealso{
  \code{\link{transport}}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
