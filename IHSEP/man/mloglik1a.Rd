\name{mloglik1a}
\alias{mloglik1a}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minus loglikelihood of an IHSEP model
}
\description{
Calculates the minus loglikelihood of an IHSEP model with given
baseline inensity function \eqn{\nu} and excitation function \eqn{g}
for event times \code{jtms} on interval \eqn{[0,TT]}.  
}
\usage{
mloglik1a(jtms, TT = max(jtms), nu, g,
          Ig = function(x) sapply(x, function(y) integrate(g, 0, y)$value),
          tol.abs = 1e-12, tol.rel = 1e-12, limit = 1000
          )

}
\arguments{
  \item{jtms}{
    A numeric vector, with values sorted in ascending order. Jump times
    to fit the inhomogeneous self-exciting point process model on. 
}
  \item{TT}{
    A scalar. The censoring time, or the terminal time for
    observation. Should be (slightly) greater than the maximum of \code{jtms}.
}
  \item{nu}{
    A (vectorized) function. The baseline intensity function. 
}
  \item{g}{
    A (vectorized) function. The excitation function. 
}
  \item{Ig}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the excitation function from 0 to \code{t}. 
  }
  \item{tol.abs}{
    A small positive number. The tolerance of the absolute error in the
    numerical integral of \eqn{\nu}. 
}
  \item{tol.rel}{
    A small positive number. The tolerance of the relative error in the
    numerical integral of \eqn{\nu}. 
}
  \item{limit}{
    An (large) positive integer. The maximum number of subintervals
    allowed in the adaptive quadrature method to find the numerical
    integral of \eqn{\nu}. 
}
  
}
\details{
  This version of the mloglik function uses external C code to speedup
  the calculations. Otherwise it is the same as the \code{mloglik0}
  function. 
}

\value{
  The value of the negative log-liklihood. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
  \code{mloglik0}
}

\examples{
## simulated data of an IHSEP on [0,1] with baseline intensity function
## nu(t)=200*(2+cos(2*pi*t)) and excitation function
## g(t)=8*exp(-16*t)
data(asep)

## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1a(tms,TT=1,nu=function(x)200*(2+cos(2*pi*x)),
          g=function(x)8*exp(-16*x),Ig=function(x)8/16*(1-exp(-16*x)))
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
\dontrun{
system.time(est <- optim(c(400,200,2*pi,8,16),
                         function(p){
                           mloglik1a(jtms=tms,TT=1,
                                     nu=function(x)p[1]+p[2]*cos(p[3]*x),
                                     g=function(x)p[4]*exp(-p[5]*x),
                                     Ig=function(x)p[4]/p[5]*(1-exp(-p[5]*x)))
                         },
                         hessian=TRUE,control=list(maxit=5000,trace=TRUE))
            )
## point estimate by MLE
est$par
## standard error estimates:
diag(solve(est$hessian))^0.5
}
}
\keyword{ point process }
\keyword{ likelihood }
