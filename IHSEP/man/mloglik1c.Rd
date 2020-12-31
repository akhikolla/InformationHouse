\name{mloglik1c}
\alias{mloglik1c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minus loglikelihood of an IHSEP model
}
\description{
Calculates the minus loglikelihood of an IHSEP model with given
baseline inensity function \eqn{\nu} and excitation function \eqn{g(x)=a
exp(-b x)} for event times \code{jtms} on interval \code{[0,TT]}.  
}
\usage{
mloglik1c(jtms, TT, nu, gcoef, Inu)
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
  \item{gcoef}{
    A numeric vector (of two elements), giving the parameters (a,b) of the
    exponential excitation function g(x)=a*exp(-b*x). 
}
  \item{Inu}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the baseline intensity function \eqn{\nu} from 0 to \code{t}. 
}

}
\details{
  This version of the mloglik function uses external C code to speedup
  the calculations. When given the analytical form of \code{Inu} or a
  quickly calculatable \code{Inu}, it should be (substantially)
  faster than \code{mloglik1a} when calculating the (minus log)
  likelihood when the excitation function is exponential. Otherwise it
  is the same as \code{mloglik0}, \code{mloglik1a}, \code{mloglik1b}. 
}

\value{
  The value of the negative log-liklihood. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
  \code{mloglik0}, \code{mloglik1a} and \code{mloglik1b}
}

\examples{
## simulated data of an IHSEP on [0,1] with baseline intensity function
## nu(t)=200*(2+cos(2*pi*t)) and excitation function
## g(t)=8*exp(-16*t)
data(asep)

## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1c(tms,TT=1,nu=function(x)200*(2+cos(2*pi*x)),
          gcoef=8*1:2,
          Inu=function(y)integrate(function(x)200*(2+cos(2*pi*x)),0,y)$value)
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
\dontrun{
system.time(est <- optim(c(400,200,2*pi,8,16),
                         function(p){
                           mloglik1c(jtms=tms,TT=1,
                                     nu=function(x)p[1]+p[2]*cos(p[3]*x),
                                     gcoef=p[-(1:3)],
                                     Inu=function(y){
                                      integrate(function(x)p[1]+p[2]*cos(p[3]*x),
                                                0,y)$value
                                     })
                         },hessian=TRUE,control=list(maxit=5000,trace=TRUE))
            )
## point estimate by MLE
est$par
## standard error estimates:
diag(solve(est$hessian))^0.5
}
}
\keyword{ point process }
\keyword{ likelihood }
