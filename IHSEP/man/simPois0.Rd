\name{simPois0}
\alias{simPois0}
\title{
Simulate a Poisson process
}
\description{
Simulate an (inhomogeneous) Poisson process with a given intensity/rate 
function over the interval \eqn{[0,T]}
}
\usage{
simPois0(int, cens = 1, int.M = optimize(int, c(0, cens), maximum = TRUE)$obj * 1.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{int}{
    A (vectorized) positive function. The intensity/rate function. 
}
  \item{cens}{
    A positive scalar. The censoring time, or time of termination of
    observations, \eqn{T}. 
}
  \item{int.M}{
    A positive scalar. Maximum value of the intensity function over \eqn{[0,T]},
    or a value larger than this. 
}
}
\details{
  The function works by first generating a Poisson process with constant
  rate \code{int.M} oever \eqn{[0,T]}, and then thinning the process
  with retention probability function
  \deqn{p(x)=\code{int}(x)/\code{int.M}}. 
}
\value{
  A numerical vector giving the event/jump times of the Poisson process
  in \eqn{[0,T]}, in ascending order. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
\code{simPois0}
}
\examples{
  aPP <- simPois(int=function(x)200*(2+cos(2*pi*x)),cens=1,int.M=600)
}
\keyword{datagen}
\keyword{point process}
\keyword{Poisson process}
