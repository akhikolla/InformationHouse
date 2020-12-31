\name{ivcoxph}

\alias{ivcoxph}

\title{
Instrumental variable estimation of the causal exposure effect in 
Cox proportional hazards (PH) models 
}

\description{
\code{ivcoxph} performs instrumental variable estimation of the causal exposure effect in 
Cox PH models with individual-level data. Below, \eqn{Z}, \eqn{X}, and 
\eqn{T} are the instrument, the exposure, and the outcome, respectively. 
\eqn{L} is a vector of covariates that we wish to control for in the analysis; 
these would typically be confounders for the instrument and the outcome.
} 

\usage{
ivcoxph(estmethod, X, fitZ.L=NULL, fitX.LZ=NULL, fitX.L=NULL, fitT.LX=NULL, 
  fitT.LZX=NULL, data, formula=~1, ctrl=FALSE, clusterid=NULL, t=NULL, 
  vcov.fit=TRUE, ...)
}

\arguments{
\item{estmethod}{
a string specifying the desired estimation method; either \code{"ts"} for two-stage
estimation, or \code{"g"} for G-estimation. 
}
  \item{X}{
a string specifying the name of the exposure \eqn{X} in \code{data}. This 
is not needed if \code{fitX.LZ} is specified. 
}
  \item{fitZ.L}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
in the \pkg{stats} package. This is a fitted GLM for \eqn{E(Z|L)}. If there are no covariates,
then \code{fitZ.L} may be specified as a model with an intercept only. This argument
is not used when \code{estmethod="ts"}. 
}
\item{fitX.LZ}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \pkg{stats} package. This is a fitted GLM for \eqn{E(X|L,Z)}. 
}
\item{fitX.L}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \pkg{stats} package. This is a fitted GLM for \eqn{E(X|L)}. If there are no covariates,
then \code{fitX.L} may be specified as a model with an intercept only. This argument
is not used when \code{estmethod="ts"}. 
}
 \item{fitT.LX}{
an object of class \code{"coxph"}, as 
  returned by the \code{coxph} function in the \pkg{survival} package. 
  This is a fitted Cox PH model for \eqn{\lambda(t|L,X)}. 
  This argument is not used when \code{estmethod="g"}.  
} 
\item{fitT.LZX}{
either an object of class \code{"coxph"} 
  or an object of class \code{"survfit"}, as returned by the \code{coxph} function 
  in the \pkg{survival} package. This is a fitted Cox PH model for 
  \eqn{\lambda(t|L,Z,X)} or a non-parametric model for \eqn{S(t|L,Z,X)},
  respectively. This argument is not used when \code{estmethod="ts"}.  
} 
  \item{data}{
a data frame containing the variables in the model. The covariates, instrument,
exposure and outcome can have arbitrary names, e.g. they don't need to 
be called \code{L}, \code{Z}, \code{X} and \code{T}.
}  
  \item{formula}{
  an object of class \code{"formula"}, with no left-hand side. This specifies
the causal interaction terms \eqn{m(L)}; see `Details'. Defaults to \code{~1}, i.e. 
main effect only. This argument is not used when \code{estmethod="ts"}.  
}
  \item{ctrl}{
logical. Should the control function \eqn{R=X-\hat{X}} be used when re-fitting 
\code{fitT.LX}? This argument is not used when \code{estmethod="g"}.   
}
  \item{clusterid}{
an optional string containing the name of a cluster identification variable when 
data are clustered. Specifying \code{clusterid} corrects the standard errors
but does not affect the estimates. 
}
 \item{t}{
 a numeric scalar specifying the time point at which to solve the estimating 
equation when \code{estmethod="g"}; see `Details'. If not specified, then the estimating equation is solved 
at the optimal value of \code{t}, defined as the value that 
minimizes \eqn{trace\{var(\hat{\psi})\}}; see Martinussen et al (2017).
This argument is not used when \code{estmethod="ts"}. 
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
\code{ivcoxph} estimates the parameter \eqn{\psi} in the causal Cox PH model
\deqn{\textrm{log}\{\lambda(t|L,Z,X)\}-\textrm{log}\{\lambda_0(t|L,Z,X)\}=m^T(L)X\psi.}
Here, \eqn{\lambda_0(t|L,Z,X)} is counterfactual hazard function,
had the exposure been set to 0. The vector function \eqn{m(L)} contains interaction terms
 between \eqn{L} and \eqn{X}. If \code{estmethod="ts"}, then these are specified 
 implicitly through the model \code{fitT.LX}. If \code{estmethod="g"}, then these     
 are specified explicitly through the \code{formula} argument.
 
If \code{estmethod="ts"}, then two-stage estimation of \eqn{\psi} is performed.
In this case, the model \code{fitX.LZ} is used to construct predictions
\eqn{\hat{X}=\hat{E}(X|L,Z)}. These predictions are subsequently used to re-fit 
the model \code{fitT.LX}, with \eqn{X} replaced with \eqn{\hat{X}}. The obtained 
coefficient(s) for \eqn{\hat{X}} in the re-fitted model is the two-stage estimator of \eqn{\psi}. 

If \code{estmethod="g"}, then G-estimation of \eqn{\psi} is performed. In this case,
the estimator is obtained as the solution to the estimating equation
\deqn{H(\psi)=\sum_{i=1}^n\hat{d}(L_i,Z_i)h_i(\psi;t)=0,}
where
\deqn{h_i(\psi;t)=\hat{S}(t|L_i,Z_i,X_i)^{\textrm{exp}\{-m^T(L_i)\psi X_i\}}.}
The estimated function \eqn{\hat{d}(L,Z)} is chosen so that the true function
has conditional mean 0, given \eqn{L}; \eqn{E\{d(L,Z)|L\}=0}. 
The specific form of \eqn{\hat{d}(L,Z)} is determined by the user-specified models. 
If \code{fitX.LZ} and \code{fitX.L} are specified, then \eqn{\hat{d}(L,Z)=m(L)\{\hat{E}(X|L,Z)-\hat{E}(X|L)\}}, 
where \eqn{\hat{E}(X|L,Z)} and \eqn{\hat{E}(X|L)} are obtained from \code{fitX.LZ}
and \code{fitX.L}, respectively. If these are not specified, then \eqn{\hat{d}(L,Z)=m(L)\{Z-\hat{E}(Z|L)\}}, 
where \eqn{\hat{E}(Z|L)} is obtained from \code{fitZ.L}, which then must be specified. 
The estimating equation is solved at the value of \eqn{t} specified by the argument \code{t}.
\eqn{\hat{S}(t|L_i,Z_i,X_i)} is an estimate of \eqn{S(t|L_i,Z_i,X_i)} obtained 
from the model \code{fitT.LZX}. 
 
}

\value{
\code{ivcoxph} returns an object of class \code{"ivcoxph"}, which inherits from 
class \code{"ivmod"}. An object of class \code{"ivcoxph"} is a list containing 
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
\code{fitT.LX}. If \code{estmethod="g"}, then the first columns correspond to \eqn{\psi}, 
and the remaining columns correspond the parameters estimated by \code{fitZ.L},
\code{fitX.LZ}, \code{fitX.L} and \code{fitT.LZX}, whichever were used in the
estimation process.
  }
\item{d.estfun}{
  the jacobian matrix of \code{colMeans(estfunall)}.    
  }
\item{converged}{
  logical. Was a solution found to the estimating equations?   
}
\item{fitT.LX}{
  the re-fitted model \code{fitT.LX} used in the estimation process when \code{estmethod="ts"}.
  This element is NULL when \code{estmethod="g"}.  
}
\item{t}{
  the value of \code{t} used in the estimation process. This element is NULL 
  when \code{estmethod="ts"}. 
}

}

\references{

Martinussen T., Sorensen D.D., Vansteelandt S. (2019). 
Instrumental variables estimation under a structural Cox model.
\emph{Biostatistics} \bold{20}(1), 65-79.

Sjolander A., Martinussen T. (2019). Instrumental variable estimation with the R package ivtools.
\emph{Epidemiologic Methods} \bold{8}(1), 1-20.

Tchetgen Tchetgen E.J., Walter S., Vansteelandt S., Martinussen T., Glymour M. (2015).
Instrumental variable estimation in a survival context. \emph{Epidemiology} \bold{26}(3), 402-410. 

}

\author{
Arvid Sjolander. 
}

\note{
\code{ivcoxph} allows for weights. However, these are defined implicitly 
through the input models. Thus, when models are used as input to \code{ivcoxph}, 
these models have to be fitted with the same weights. When \code{estmethod="g"}
the weights are taken from \code{fitX.LZ}, if specified by the user. If \code{fitX.LZ} is not
specified then the weights are taken from \code{fitZ.L}. Hence, if weights are used,
then either \code{fitX.LZ} or \code{fitZ.L} must be specified.   
}

\examples{

require(survival)

set.seed(9)

##Note: the parameter values in the examples below are chosen to make 
##Y0 independent of Z, which is necessary for Z to be a valid instrument.

n <- 10000
psi0 <- 0.5
Z <- rbinom(n, 1, 0.5)
X <- rbinom(n, 1, 0.7*Z+0.2*(1-Z))
m0 <- exp(0.8*X-0.41*Z) #T0 independent of Z at t=1
T <- rexp(n, rate=exp(psi0*X+log(m0)))
C <- rexp(n, rate=exp(psi0*X+log(m0))) #50% censoring
d <- as.numeric(T<C)
T <- pmin(T, C)
data <- data.frame(Z, X, T, d)

#two-stage estimation
fitX.LZ <- glm(formula=X~Z, data=data)
fitT.LX <- coxph(formula=Surv(T, d)~X, data=data)
fitIV <- ivcoxph(estmethod="ts", fitX.LZ=fitX.LZ, fitT.LX=fitT.LX, data=data, 
  ctrl=TRUE)
summary(fitIV)

#G-estimation with non-parametric model for S(t|L,Z,X) and model for Z
fitZ.L <- glm(formula=Z~1, data=data)
fitT.LZX <- survfit(formula=Surv(T, d)~X+Z, data=data)
fitIV <- ivcoxph(estmethod="g", X="X", fitZ.L=fitZ.L, fitT.LZX=fitT.LZX, 
  data=data, t=1)
summary(fitIV)

#G-estimation with Cox model for \lambda(t|L,Z,X) and model for Z
fitZ.L <- glm(formula=Z~1, data=data)
fitT.LZX <- coxph(formula=Surv(T, d)~X+X+X*Z, data=data)
fitIV <- ivcoxph(estmethod="g", X="X", fitZ.L=fitZ.L, fitT.LZX=fitT.LZX, 
  data=data, t=1)
summary(fitIV)

}
