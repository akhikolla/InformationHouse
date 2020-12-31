\name{semidiscrete}
\alias{semidiscrete}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find Optimal Transport Partition Between pgrid and wpp.
}
\description{
  Given an object \code{a} of class \code{\link[=pgrid-object]{pgrid}} specifying an image and an object \code{b}
  of class \code{\link[=wpp-object]{wpp}} specifiying a more flexible mass distribution at finitely many points,
  find the partition of the image (and hence the optimal transport map) that minimizes the total transport cost
  for going from \code{a} to \code{b}.
}
\usage{
  semidiscrete(a, b, p = 2, method = c("aha"), control = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
  	an object of class \code{\link[=pgrid-object]{pgrid}} usually representing an image or the discretization of a measure.
}
  \item{b}{
    an object of class \code{\link[=wpp-object]{wpp}} usually having the same total mass as \code{a}.   	
}
  \item{p}{
    the power \eqn{\geq 1}{>=1} to which the Euclidean distance between points is taken in order
    to compute costs. Only \eqn{p \in \{1,2\}}{p in {1,2}} is implemented.
}
  \item{method}{
    the name of the algorithm to use. Currently only \code{aha} is supported.
}
  \item{control}{
  	a named list of parameters for the chosen method or the result of a call to \code{\link{trcontrol}}. Currently only
  	the parameters \code{factr} and \code{maxit} can be set.
  }
  \item{...}{
  	currently without effect.
  }
}
\details{
  This is a wrapper for the functions \code{\link{aha}} and \code{\link{semidiscrete1}}. In the former the
  Aurenhammer--Hoffmann--Aronov (1998) method for \eqn{p=2} is implemented in the multiscale variant presented
  in Mérigot (2011). In the latter an adapted Aurenhammer--Hoffmann--Aronov method for \eqn{p=1} is used that
  was presented in Hartmann and Schuhmacher (2018).
  
  The present function is automatically called by \code{\link{transport}} if the
  first argument is of class \code{pgrid} and the second argument is of class \code{wpp}.
}
\value{
  An object describing the optimal transport partition for \code{a} and \code{b}.
  
  If \code{p=1} an object of class \code{apollonius_diagram} having components \code{sites} and \code{weights},
  as well as (optionally) \code{wasserstein_dist} and \code{ret_code} (the return code from the call to
  \code{\link{semidiscrete1}}).
  
  If \code{p=2} an objectof class \code{power_diagram} having components \code{sites} and \code{cells},
  as well as (optionally) \code{wasserstein_dist}. \code{sites} is here a data.frame with columns \code{xi},
  \code{eta} and \code{w} (the weights for the power diagram). \code{cells} is a list with as many
  2-column matrix components as there are sites, each describing the \eqn{x}- and \eqn{y}-coordinates
  of the polygonal cell associated with the corresponding site or \code{NULL} if the cell of the site is empty.
  
  Plotting methods exist for objects of class \code{apollonius_diagram}, \code{power_diagram} and
  for \link[=plot.pgrid]{optimal transport maps represented by either of the two}.
}
\references{
F. Aurenhammer, F. Hoffmann and B. Aronov (1998). Minkowski-type theorems and least-squares clustering. Algorithmica 20(1), 61--76.

V. Hartmann and D. Schuhmacher (2017). Semi-discrete optimal transport --- the case p=1. Preprint \href{https://arxiv.org/abs/1706.07650}{arXiv:1706.07650}

M. Karavelas and M. Yvinec. 2D Apollonius Graphs 
(Delaunay Graphs of Disks). In CGAL User and Reference Manual.
CGAL Editorial Board, 4.12 edition, 2018
            
Q. Mérigot (2011). A multiscale approach to optimal transport. Computer Graphics Forum 30(5), 1583--1592. \href{http://dx.doi.org/10.1111/j.1467-8659.2011.02032.x}{doi:10.1111/j.1467-8659.2011.02032.x}

Naoaki Okazaki (2010). libLBFGS: a library of Limited-memory
Broyden-Fletcher-Goldfarb-Shanno (L-BFGS). Version 1.10
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de} \cr  
  Björn Bähre \email{bjobae@gmail.com} \cr  
  Valentin Hartmann \email{valentin.hartmann@epfl.ch}
}
\note{
For \code{p=1} this function requires the Computational Geometry Algorithms Library (CGAL), available at \url{https://www.cgal.org}. Adapt the file src/Makevars according to the instructions given there and re-install from source.

Internally the code from liblbfgs 1.10 by Naoaki Okazaki (2010) is used.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=plot.pgrid]{plot}}, \code{\link{transport}}, \code{\link{aha}}, \code{\link{semidiscrete1}}
}
\examples{
##  See examples for function transport}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{semi-discrete}
\keyword{transport partition}
\keyword{transport tesselation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
