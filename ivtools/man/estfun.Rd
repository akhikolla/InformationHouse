\name{estfun}

\alias{estfun}

\title{
Computes the estimating function sum for \code{"ivmod"}  
objects, fitted with \code{estmethod="g"}.
}

\description{
\code{estfun} computes the estimating function \eqn{H(\psi)}
for a \code{"ivmod"} object, fitted with \code{estmethod="g"}, 
for a range of values of \eqn{\psi}. The \code{estfun} is not implemented 
for \code{"ivah"} objects, since G-estimation in additive hazards models
is based on a recursive estimation technique, and not standard estimating equations.
}     

\usage{
estfun(object, lower, upper, step)
}

\arguments{
\item{object}{
an object of class \code{"ivmod"}, fitted with \code{estmethod="g"}.   
}
  \item{lower}{
an optional vector of lower values for \eqn{\psi}. Defaults to \eqn{\psi-0.5}. 
}
\item{upper}{
an optional vector of upper values for \eqn{\psi}. Defaults to \eqn{\psi+0.5}. 
}
\item{step}{
an optional vector of steps between \code{lower} and \code{upper}. Defaults to 0.01
for each element of \eqn{\psi}. 
}

}

\details{ 
\code{estfun} may be useful for visual inspection of the estimating 
function, to make sure that a solution to the estimating equation 
\deqn{H(\psi)=0}
was found, see `Examples'. For the \eqn{i}:th element of \eqn{\psi}, the estimating
function sum is computed for a range of values within (\code{lower[i]}, \code{upper[i]}), at the 
G-estimate of the remaining elements of \eqn{\psi}.  
 
}

\value{

An object of class \code{"estfun"} is a list containing 
\item{f}{
  a named list of matricies; one matrix for each element of \eqn{\psi}. 
  The first column of the \eqn{i}:th matrix contains the values for the \eqn{i}:th
  element of \eqn{\psi} at which the estimating function sum is computed, the 
  second column contains the values of the estimating function sum. 
  }

\item{est}{
  the G-estimate of \eqn{\psi}.
}
}


\references{
Burgess S, Granell R, Palmer TM, Sterne JA, Didelez V. (2014). 
Lack of identification in semiparametric instrumental variable models with 
binary outcomes. \emph{American Journal of Epidemiology} \bold{180}(1), 111-119.
 
Vansteelandt S., Bowden J., Babanezhad M., Goetghebeur E. (2011). 
On instrumental variables estimation of causal odds ratios.
\emph{Statistical Science} \bold{26}(3), 403-422.

}


\author{
Arvid Sjolander. 
}

\examples{

set.seed(9)

##Note: the parameter values in the examples below are chosen to make 
##Y0 independent of Z, which is necessary for Z to be a valid instrument.

n <- 1000
psi0 <- 0.5
psi1 <- 0.2

##---Example 1: linear model and interaction between X and L---

L <- rnorm(n)
Z <- rnorm(n, mean=L)
X <- rnorm(n, mean=Z)
m0 <- X-Z+L 
Y <- rnorm(n, mean=psi0*X+psi1*X*L+m0)
data <- data.frame(L, Z, X, Y)

#G-estimation
fitZ.L <- glm(formula=Z~L, data=data)
fitIV <- ivglm(estmethod="g", X="X", Y="Y", fitZ.L=fitZ.L, data=data, 
  formula=~L, link="identity")
summary(fitIV)
H <- estfun(fitIV)
plot(H)

##---Example 2: logistic model and no covariates--- 

Z <- rbinom(n, 1, 0.5)
X <- rbinom(n, 1, 0.7*Z+0.2*(1-Z)) 
m0 <- plogis(1+0.8*X-0.39*Z)
Y <- rbinom(n, 1, plogis(psi0*X+log(m0/(1-m0)))) 
data <- data.frame(Z, X, Y)

#G-estimation
fitZ.L <- glm(formula=Z~1, data=data)
fitY.LZX <- glm(formula=Y~X+Z+X*Z, family="binomial", data=data)
fitIV <- ivglm(estmethod="g", X="X", fitZ.L=fitZ.L, fitY.LZX=fitY.LZX, 
  data=data, link="logit")
summary(fitIV)
H <- estfun(fitIV)
plot(H)

}
