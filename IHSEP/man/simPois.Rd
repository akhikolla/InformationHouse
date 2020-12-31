\name{simPois}
\alias{simPois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a Poisson process
}
\description{
Simulate an (inhomogeneous) Poisson process with a given intensity/rate 
function over the interval \eqn{[0,T]}
}
\usage{
simPois(int=NULL,cens=1,int.M=optimize(int,c(0,cens),maximum=TRUE)$obj*1.1,
        B=max(as.integer(sqrt(int.M * cens)),10), exp.int=FALSE,par=c(1,2))
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
    or a larger value. 
}
  \item{B}{
    An integer. The block size to be used in generating exponential
    random variables in blocks. 
  }
  \item{exp.int}{
    Logical. Set to TRUE, if the intensity function is exponential,
    i.e. a*exp(-b*x). If set to TRUE, the parameters a and b should also
    be supplied via the argument par.
  }
  \item{par}{
    A numerical vector of two elements, giving the parameters a and b of
    the exponential intensity function. The values are not ignored if
    exp.int is set to FALSE. 
  }
}
\details{
  The function works by first generating a Poisson process with constant
  rate \code{int.M} oever \eqn{[0,T]}, and then thinning the process
  with retention probability function
  \deqn{p(x)=\code{int}(x)/\code{int.M}}. 

  When generating the homoneous Poisson process, it first generates
  about \eqn{\Lambda+1.96*\sqrt{\Lambda}} exponential variables, then,
  if the exponential variables are not enough (their sum has not reached
  the censoring time \eqn{T} yet), generates exponential variables in
  blocks of size \code{B} until the total of all the generated
  exponentials exceeds \eqn{T}. Then \code{cumsum}s of the exponentials
  that are smaller than or equal to \eqn{T} are retained as the event
  times of the homogeneous Poisson process. This method apparantly does
  not produce tied event times. 
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (int, cens = 1, int.M = optimize(int, c(0, cens), maximum = TRUE)$obj * 
    1.1, B = max(as.integer(sqrt(int.M * cens)), 10)) 
{
    tms <- rexp(as.integer(int.M * cens + 2 * sqrt(int.M * cens)), 
        rate = int.M)
    while (sum(tms) < cens) tms <- c(tms, rexp(B, rate = int.M))
    cumtms <- cumsum(tms)
    tms <- cumtms[cumtms <= cens]
    N <- length(tms)
    if (N == 0) 
        return(numeric(0))
    tms[as.logical(mapply(rbinom, n = 1, size = 1, prob = int(tms)/int.M))]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
