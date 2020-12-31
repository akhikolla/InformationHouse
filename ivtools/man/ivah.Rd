\name{ivah}

\alias{ivah}

\title{
Instrumental variable estimation of the causal exposure effect in 
additive hazards (AH) models 
}

\description{
\code{ivah} performs instrumental variable estimation of the causal exposure effect in 
AH models with individual-level data. Below, \eqn{Z}, \eqn{X}, and 
\eqn{T} are the instrument, the exposure, and the outcome, respectively. 
\eqn{L} is a vector of covariates that we wish to control for in the analysis; 
these would typically be confounders for the instrument and the outcome.
} 

\usage{
ivah(estmethod, X, T, fitZ.L=NULL, fitX.LZ=NULL, fitT.LX=NULL, data,  
  ctrl=FALSE, clusterid=NULL, event, max.time, max.time.psi, n.sim=100, 
  vcov.fit=TRUE, ...)
}

\arguments{
\item{estmethod}{
a string specifying the desired estimation method; either \code{"ts"} for two-stage
estimation, or \code{"g"} for G-estimation. 
}
  \item{X}{
a string specifying the name of the exposure \eqn{X} in \code{data}. This is not needed if 
\code{fitX.LZ} is specified.
}
  \item{T}{
a string specifying the name of the follow-up time \eqn{T} in \code{data}. This is not needed if 
\code{fitT.LZ} is specified.
}
  \item{fitZ.L}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
in the \pkg{stats} package. This is a fitted GLM for \eqn{E(Z|L)}. If there are no covariates,
then \code{fitZ.L} may be specified as a model with an intercept only. This argument
is not used when \code{estmethod="ts"}. 
}
\item{fitX.LZ}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \pkg{stats} package. This is a fitted GLM for \eqn{E(X|L,Z)}. This argument
is not used when \code{estmethod="g"}. 
}
 \item{fitT.LX}{
If \code{estmethod="ts"}, then this is an object of class \code{"ah"}, as 
  returned by the \code{ah} function in the \pkg{ivtools} package. 
  In this case it is a fitted AH model for \eqn{\lambda(t|L,X)}. 
  This argument is not used when \code{estmethod="g"}.  
} 
  \item{data}{
a data frame containing the variables in the model. The covariates, instrument,
exposure and outcome can have arbitrary names, e.g. they don't need to 
be called \code{L}, \code{Z}, \code{X} and \code{T}.
}  
  \item{ctrl}{
logical. Should the control function \eqn{R=X-\hat{X}} be used when re-fitting 
\code{fitY}? This argument is not used when \code{estmethod="g"}.   
}
  \item{clusterid}{
an optional string containing the name of a cluster identification variable when 
data are clustered. Specifying \code{clusterid} corrects the standard errors
but does not affect the estimates. This argument is not used when 
\code{estmethod="g"}, since correction for clustered data is currently not implemented
for G-estimation. 
}
  \item{event}{
a string specifying the name of the status indicator, 0="no event", 1="event".
This argument is not used when \code{estmethod="ts"}. 
}

  \item{max.time}{
optional follow-up for estimating \eqn{B(t)} with G-estimation. 
Defaults to maximal observed follow-up time in \code{data}.
This argument is not used when \code{estmethod="ts"}. 
}
  \item{max.time.psi}{
optional follow-up for estimating \eqn{\psi} with G-estimation. 
Defaults to maximal observed follow-up time in \code{data}.
This argument is not used when \code{estmethod="ts"}.
}
  \item{n.sim}{
optional number of resamplings for testing goodness-of-fit of constant effects model
for G-estimation. Defaults to 100. This argument is not used when \code{estmethod="ts"}.  
}
 \item{vcov.fit}{
 logical. Should the variance-covariance matrix be computed? 
}
  \item{...}{
  optional arguments passed on to the \code{nleqslv} function, which is used to 
  solve the estimating equations when \code{estmethod="g"}. See the help pages 
  for \code{nleqslv}. This argument is not used when \code{estmethod="ts"}. 
  } 

}

\details{ 
The \code{ivah} estimates different parameters, depending on whether
\code{estmethod="ts"} or \code{estmethod="g"}. If \code{estmethod="ts"}, then 
\code{ivah} uses two-stage estimation to estimate the parameter \eqn{\psi} in the causal AH model
\deqn{\lambda(t|L,Z,X)-\lambda_0(t|L,Z,X)=m^T(L)X\psi.}
Here, \eqn{\lambda_0(t|L,Z,X)} is counterfactual hazard function,
had the exposure been set to 0. The vector function \eqn{m(L)} contains interaction terms
 between \eqn{L} and \eqn{X}. These are specified 
 implicitly through the model \code{fitY}. The model \code{fitX.LZ} is used to 
 construct predictions \eqn{\hat{X}=\hat{E}(X|L,Z)}. These predictions are 
 subsequently used to re-fit the model \code{fitY}, with \eqn{X} replaced with 
 \eqn{\hat{X}}. The obtained coefficient(s) for \eqn{X} is the two-stage 
 estimator of \eqn{\psi}.  
 
If \code{estmethod="g"}, then \code{ivah} uses G-estimation to estimate the function 
\eqn{B(t)} in the causal AH model
\deqn{\lambda(t|L,Z,X)-\lambda_0(t|L,Z,X)=XdB(t).} 
It also delivers an estimate of \eqn{dB(t)} assuming that this function is 
constant across time (=\eqn{\psi}).
}

\value{
\code{ivah} returns an object of class \code{"ivah"}, which inherits from 
class \code{"ivmod"}. An object of class \code{"ivah"} is a list containing 
\item{call}{
  the matched call. 
  }
\item{input}{
  \code{input} is a list containing all input arguments     
  }
\item{est}{
  a vector containing the estimate of \eqn{\psi}.  
  }
\item{vcov}{
  the variance-covariance matrix for the estimate of \eqn{\psi}, 
obtained with the sandwich formula. 
  }
\item{estfunall}{
  a matrix of all subject-specific contributions to the estimating functions used in the estimation process.
One row for each subject, one column for each parameter. If \code{estmethod="ts"},
then the first columns correspond to the parameters estimated by \code{fitX.LZ}, and
the last columns correspond to the parameters estimated by the re-fitted model 
\code{fitY}. If \code{estmethod="g"}, then \code{estfunall} is NULL. 
  }
\item{d.estfun}{
  the jacobian matrix of \code{colMeans(estfun)}. If \code{estmethod="g"}, 
    then \code{d.estfun} is NULL.    
  }
\item{converged}{
  logical. Was a solution found to the estimating equations?   
}
\item{fitY}{
  the re-fitted model \code{fitY} used in the estimation process when \code{estmethod="ts"}.
  This element is NULL when \code{estmethod="g"}.  
}
\item{stime}{
  the ordered event times within (0,max.time). This element is NULL when \code{estmethod="ts"}. 
}
\item{B}{
  the estimate of \eqn{B(t)}. This element is NULL when \code{estmethod="ts"}.   
}
\item{se_B}{
  the standard error of the estimate of \eqn{B(t)}. This element is NULL when \code{estmethod="ts"}.
}
\item{pval_0}{
  p-value corresponding to supremum test of the null \eqn{B(t)=0}. 
  This element is NULL when \code{estmethod="ts"}.
}
\item{eps_B}{
  the iid-decomposition of \eqn{\sqrt{n}(\hat{B}(t) - B(t))}. 
  This element is NULL when \code{estmethod="ts"}.
}
\item{pval_psi}{
  p-value corresponding to the null \eqn{\psi=0}. 
  This element is NULL when \code{estmethod="ts"}.
}
\item{pval_GOF_sup}{
  p-value corresponding to supremum test of the null \eqn{B(t)=\psi}.  
  This element is NULL when \code{estmethod="ts"}.
}
\item{pval_GOF_CvM}{
  as pval_GOF_sup but now based on the Cramer Von Mises test statistic.  
  This element is NULL when \code{estmethod="ts"}.
}
\item{GOF.resamp}{
  a matrix with first row the ordered jump times in (0,max.time.bet),
 second row the observed test process, and the remaining rows are 
 50 processes sampled under the null. 
  This element is NULL when \code{estmethod="ts"}.
}
}

\references{

Martinussen T., Vansteelandt S., Tchetgen Tchetgen E.J., Zucker D.M. (2017). 
Instrumental variables estimation of exposure effects on a time-to-event 
endpoint using structural cumulative survival models. 
\emph{Epidemiology} \bold{73}(4): 1140-1149.

Sjolander A., Martinussen T. (2019). Instrumental variable estimation with the R package ivtools.
\emph{Epidemiologic Methods} \bold{8}(1), 1-20.

Tchetgen Tchetgen E.J., Walter S., Vansteelandt S., Martinussen T., Glymour M. (2015).
Instrumental variable estimation in a survival context. \emph{Epidemiology} \bold{26}(3): 402-410. 

}

\author{
Arvid Sjolander and Torben Martinussen. 
}

\note{
\code{ivah} allows for weights. However, these are defined implicitly 
through the input models. Thus, when models are used as input to \code{ivah}, 
these models have to be fitted with the same weights.

Left-truncation and correction of standard errors for clustered data are currently not
implemented when \code{estmethod="g"}.     
}

\examples{

require(ahaz)

set.seed(9)

n <- 1000
psi0 <- 0.2
psi1 <- 0.0

U <- runif(n)
L <- runif(n)
Z <- rbinom(n, 1, plogis(-0.5+L)) 
X <- runif(n, min=Z+U, max=2+Z+U)
T <- rexp(n, rate=psi0*X+psi1*X*L+0.2*U+0.2*L)
C <- 5  #administrative censoring at t=5
d <- as.numeric(T<C)
T <- pmin(T, C) 
data <- data.frame(L, Z, X, T, d)
#break ties
data$T <- data$T+rnorm(n=nrow(data), sd=0.001)

#two-stage estimation
fitX.LZ <- glm(formula=X~Z+L, data=data)
fitT.LX <- ah(formula=Surv(T, d)~X+L+X*L, data=data)
fitIV <- ivah(estmethod="ts", fitX.LZ=fitX.LZ, fitT.LX=fitT.LX, data=data, 
  ctrl=TRUE)
summary(fitIV)

#G-estimation
fitZ.L <- glm(formula=Z~L, family="binomial", data=data)
fitIV <- ivah(estmethod="g", X="X", T="T", fitZ.L=fitZ.L, data=data,  
  event="d", max.time=4, max.time.psi=4, n.sim=100)
summary(fitIV)
plot(fitIV)



}
