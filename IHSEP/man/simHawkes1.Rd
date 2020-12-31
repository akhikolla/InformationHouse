\name{simHawkes1}
\alias{simHawkes1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a Hawkes process, or Self-exciting point process
}
\description{
Simulate an (inhomogeneous) self-exciting process with given background
intensity and excitation/fertility function. 
}
\usage{
simHawkes1(nu=NULL, g=NULL, cens = 1,
           nuM=max(optimize(nu,c(0,cens),maximum=TRUE)$obj, nu(0), nu(cens))*1.1,
           gM=max(optimize(g,c(0,cens),maximum = TRUE)$obj, g(0),
g(cens))*1.1,
           exp.g=FALSE,gp=c(1,2))
}
%- maybe also 'usage' for other objects documented here.
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
  \item{exp.g}{
    A logical. Whether the excitation function \code{g} should be treated as an
    exponential function.
  }
  \item{gp}{
    A vector of two elements, giving the two parameters a and b in the
    exponential excitation function \eqn{g(x)=a*exp(-b*x)}, which is
    used when \code{exp.g} is set to TRUE, and is ignored otherwise. 
  }
}
\details{
  The function works by simulating the birth times generation by
  generation according to inhomegenous Poisson processes with
  appropriate intensity functions (\eqn{\nu} or \eqn{g}). 
}
\value{
  A list of vectors of arrival/birth times of individuals/events of
  generations 0, 1, ... The length of the list is the total number of
  generations.   
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{simHawkes0}
}
\examples{
  asepp <- simHawkes1(nu=function(x)200*(2+cos(2*pi*x)),nuM=600,
                      g=function(x)8*exp(-16*x),gM=8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{point process}
\keyword{self-exciting}
