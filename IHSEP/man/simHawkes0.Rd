\name{simHawkes0}
\alias{simHawkes0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a Hawkes process, or Self-exciting point process
}
\description{
Simulate an (inhomogeneous) self-exciting process with given background
intensity and excitation/fertility function. 
}
\usage{
simHawkes0(nu, g, cens = 1,
           nuM=max(optimize(nu,c(0,cens),maximum=TRUE)$obj,nu(0),nu(cens))*1.1,
           gM=max(optimize(g,c(0,cens),maximum=TRUE)$obj, g(0),g(cens))*1.1)
}
\arguments{
  \item{nu}{
        A (vectorized) function. The baseline intensity function. 
}
  \item{g}{
    A (vectorized) function. The excitation function. 
}
  \item{cens}{
    A scalar. The censoring time, or the time of termination of observations. 
}
  \item{nuM}{
    A scalar. The maximum time of the baseline intensity from 0 to \code{cens}.
}
  \item{gM}{
    A scalar. The maximum time of the excitation function from 0 to \code{cens}.
}
}
\details{
  The function works by simulating the birth times generation by
  generation according to inhomegenous Poisson processes with
  appropriate intensity functions (\eqn{\nu} or \eqn{g}). 
}
\value{
  A list of vectors of arrival/birth times of individuals/events of
  generations 0, 1, .... 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\note{
  Same algorithm as in \code{simHawkes1}, though the latter might be
  more succinct and (very slightly) faster. 
}
\seealso{
  \code{simHawkes1}
}
\examples{
asepp <- simHawkes0(nu=function(x)200*(2+cos(2*pi*x)),nuM=600,
                               g=function(x)8*exp(-16*x),gM=8)
}
\keyword{datagen}
\keyword{point process}
\keyword{self-exciting}
