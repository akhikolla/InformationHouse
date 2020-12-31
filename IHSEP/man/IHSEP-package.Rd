\name{IHSEP-package}
\alias{IHSEP-package}
\alias{IHSEP}
\docType{package}
\title{
  Inhomogeneous Self-exciting Process
}
\description{
  Simulates the (inhomogeneous) Self-exciting point process (SEPP), or
  Hawkes process, on \eqn{[0,T]} with a given (possibly time-varying)
  baseline/background intensity function \eqn{nu} and excitation
  function (fertility rate function) \eqn{g}.  Or calculate the
  likelihood of an SEPP with baseline intensity \eqn{nu} and excitation
  function \eqn{g} for a given set of event times on \eqn{[0,T]}.
}
\details{
\tabular{ll}{
Package: \tab IHSEP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-12\cr
License: \tab GPL(>=2)\cr
}
}
\author{
Feng Chen <feng.chen@unsw.edu.au>

Maintainer: Feng Chen <feng.chen@unsw.edu.au>
}
\references{
Feng Chen and Peter Hall (2013). Inference for a nonstationary
self-exciting point process with an application in ultra-high frequency
financial data modeling.  Journal of Applied Probability
50(4):1006-1024.

Feng Chen and Peter Hall (2016). Nonparametric Estimation for
Self-Exciting Point Processes -- A Parsimonious Approach. Journal of
Computational and Graphical Statistics 25(1): 209-224. 

Alan G Hawkes (1971). Spectra of some self-exciting and mutually
exciting point processes. Biometrika 58(1):83-90.
}
\keyword{ package }
\examples{
\dontrun{
## simulate the times of a Poisson process on [0,1] with intensity
## function nu(t)=100*(2+cos(2*pi*t)).
tms <- simPois(int=function(x)100*(2+cos(2*pi*x)),int.M=300)
## calculate a nonparametric estimate of the intensity function
int <- lpint::lpint(tms,Tau=1)
matplot(int$x,int$y+qnorm(0.975)*outer(int$se,c(-1,0,1)),type="l",lty=c(2,1,2),
        col=1,xlab="t",ylab="nu(t)")
curve(100*(2+cos(2*pi*x)),add=TRUE,col="red")

## simulate an IHSEP on [0,1] with baseline intensity function
## nu(t)=100*(2+cos(2*pi*t)) and excitation function
## g(t)=0.5*8*exp(-8*t)
asep <- simHawkes1(nu=function(x)200*(2+cos(2*pi*x)),nuM=600,
                              g=function(x)8*exp(-16*x),gM=8)
## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1a(tms,TT=1,nu=function(x)200*(2+cos(2*pi*x)),
          g=function(x)8*exp(-16*x),Ig=function(x)8/16*(1-exp(-18*x)))
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
est <- optim(c(400,200,2*pi,8,16),
             function(p){
               mloglik1a(jtms=tms,TT=1,
                         nu=function(x)p[1]+p[2]*cos(p[3]*x),
                         g=function(x)p[4]*exp(-p[5]*x),
                         Ig=function(x)p[4]/p[5]*(1-exp(-p[5]*x)))
             },
             hessian=TRUE,control=list(maxit=5000,trace=TRUE))
## point estimate by MLE
est$par
## standard error estimates:
diag(solve(est$hessian))^0.5
}
}