\name{transport-package}
\alias{transport-package}
\docType{package}
\encoding{UTF-8}
\title{
Optimal Transport in Various Forms
}
\description{
Solve optimal transport problems. Compute Wasserstein distances (a.k.a. Kantorovitch, Fortet--Mourier, Mallows, Earth Mover's, or minimal \eqn{L_p} distances), return the corresponding transport plans, and display them graphically. Objects that can be compared include grey-scale images, (weighted) point patterns, and mass vectors.
}
\details{
\tabular{ll}{
Package: \tab transport\cr
Type: \tab Package\cr
Version: \tab 0.12-1\cr
Date: \tab 2019-08-07\cr
License: \tab GPL (>=2) \cr
LazyData: \tab yes\cr
}

The main end-user function is \code{\link{transport}}. It computes optimal transport plans between images (class \code{\link{pgrid}}), point patterns (class \code{\link{pp}}), weighted point patterns (class \code{\link{wpp}}) and mass vectors, based on various algorithms. These transport plans can be \code{\link[=plot.pgrid]{plot}}ed. The function \code{\link{wasserstein}} allows for the numerical computation of \eqn{p}-th order Wasserstein distances.

Most functions in this package are designed for data in two and higher dimensions. A quick tool for computing the \eqn{p}-th order Wasserstein distance between univariate samples is \code{\link{wasserstein1d}}.

}
\author{
  Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
  
  Björn Bähre \email{bjobae@gmail.com} (code for \code{\link{aha}}-method)
  
  Nicolas Bonneel \email{nicolas.bonneel@liris.cnrs.fr} \cr (adaptation of LEMON code for fast \code{\link[=transport]{networkflow}} method)
  
  Carsten Gottschlich \email{gottschlich@math.uni-goettingen.de} \cr (original java code for \code{\link[=transport]{shortlist}} and \code{\link[=transport]{revsimplex}} methods)
  
  Valentin Hartmann \email{valentin.hartmann@epfl.ch} (code for \code{\link[=semidiscrete1]{aha}} method for \code{p=1})
  
  Florian Heinemann \email{florian.heinemann@uni-goettingen.de} \cr (integration of \code{\link[=transport]{networkflow}} method)
  
  Bernhard Schmitzer \email{schmitzer@uni-muenster.de} (\code{\link{shielding}} method)
  
  Jörn Schrieber \email{joern.schrieber-1@mathematik.uni-goettingen.de} (\code{\link[=subwasserstein]{subsampling}} method)

Maintainer: Dominic Schuhmacher \email{dominic.schuhmacher@mathematik.uni-goettingen.de}
}
\references{
  See help page for the function \code{\link{transport}}.
}
\keyword{package}
\examples{
##  See examples for function transport
}
