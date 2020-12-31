\name{mloglik1b}
\alias{mloglik1b}
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
mloglik1b(jtms, TT = max(jtms), nu, g,
          Ig=function(x)sapply(x,function(y)integrate(g,0,y,
               rel.tol=1e-12,abs.tol=1e-12,subdivisions=1000)$value),
          Inu=function(x)sapply(x,function(y)integrate(nu,0,y)$value))
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
  \item{Inu}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the baseline intensity function \eqn{\nu} from 0 to \code{t}. 
}

}
\details{
  This version of the mloglik function uses external C code to speedup
  the calculations. When given the analytical form of \code{Inu} or a
  quickly calculatable \code{Inu}, it should be (probably slightly)
  faster than \code{mloglik1a}. Otherwise it is the same as
  \code{mloglik0} and \code{mloglik1a}. 
}

\value{
  The value of the negative log-liklihood. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
  \code{mloglik0} and \code{mloglik1a}
}

\examples{
## simulated data of an IHSEP on [0,1] with baseline intensity function
## nu(t)=200*(2+cos(2*pi*t)) and excitation function
## g(t)=8*exp(-16*t)
data(asep)

## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1b(tms,TT=1,nu=function(x)200*(2+cos(2*pi*x)),
          g=function(x)8*exp(-16*x),Ig=function(x)8/16*(1-exp(-16*x)))
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
\dontrun{
system.time(est <- optim(c(400,200,2*pi,8,16),
                         function(p){
                           mloglik1b(jtms=tms,TT=1,
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
