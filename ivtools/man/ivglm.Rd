\name{ivglm}

\alias{ivglm}

\title{
Instrumental variable estimation of the causal exposure effect in 
generalized linear models 
}

\description{
\code{ivglm} performs instrumental variable estimation of the causal exposure effect in 
generalized linear models with individual-level data. Below, \eqn{Z}, \eqn{X}, and 
\eqn{Y} are the instrument, the exposure, and the outcome, respectively. 
\eqn{L} is a vector of covariates that we wish to control for in the analysis; 
these would typically be confounders for the instrument and the outcome.
}     

\usage{
ivglm(estmethod, X, Y, fitZ.L=NULL, fitX.LZ=NULL, fitX.L=NULL, fitY.LX=NULL, 
  fitY.LZX=NULL, data, formula=~1, ctrl=FALSE, clusterid=NULL, link, vcov.fit=TRUE, 
  ...)
}

\arguments{
\item{estmethod}{
a string specifying the desired estimation method; either \code{"ts"} for two-stage
estimation, or \code{"g"} for G-estimation. 
}
  \item{X}{
a string specifying the name of the exposure \eqn{X} in \code{data}. This is not needed 
if \code{fitX.LZ} is specified. 
}
  \item{Y}{
a string specifying the name of the outcome \eqn{Y} in \code{data}. This 
is not needed if \code{fitY.LX} or \code{fitY.LZX} is specified. 
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
 \item{fitY.LX}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \pkg{stats} package. This is a fitted GLM for \eqn{E(Y|L,X)}. 
This argument is not used when \code{estmethod="g"}.  
} 
 \item{fitY.LZX}{
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \pkg{stats} package. This is a fitted GLM for \eqn{E(Y|L,Z,X)}. 
This argument is not used when \code{estmethod="ts"}. It is also not used 
when \code{estmethod="g"} and \code{link="identity"} or \code{link="log"}.  
} 
  \item{data}{
a data frame containing the variables in the model. The covariates, instrument,
exposure and outcome can have arbitrary names, e.g. they don't need to 
be called \code{L}, \code{Z}, \code{X} and \code{Y}.
}  
  \item{formula}{
  an object of class \code{"formula"}, with no left-hand side. This specifies
the causal interaction terms \eqn{m(L)}; see `Details'. Defaults to \code{~1}, i.e. 
main effect only. This argument is not used when \code{estmethod="ts"}.  
}
  \item{ctrl}{
logical. Should the control function \eqn{R=X-\hat{X}} be used when re-fitting 
\code{fitY.LX}? This argument is not used when \code{estmethod="g"}.   
}
  \item{clusterid}{
an optional string containing the name of a cluster identification variable when 
data are clustered. Specifying \code{clusterid} corrects the standard errors
but does not affect the estimates. 
}
 \item{link}{
 a string specifying the link function for the causal generalized linear model;
 see `Details'. Either \code{"identity"}, \code{"log"}, or \code{"logit"}. 
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
\code{ivglm} estimates the parameter \eqn{\psi} in the causal generalized linear model
\deqn{\eta\{E(Y|L,Z,X)\}-\eta\{E(Y_0|L,Z,X)\}=m^T(L)X\psi.}
Here, \eqn{E(Y_0|L,Z,X)} is counterfactual mean of the outcome, 
had the exposure been set to 0. The link function
 \eqn{\eta} is either the identity, log or logit link, as specified by 
 the \code{link} argument. The vector function \eqn{m(L)} contains interaction terms
 between \eqn{L} and \eqn{X}. If \code{estmethod="ts"}, then these are specified 
 implicitly through the model \code{fitY.LX}. If \code{estmethod="g"}, then these     
 are specified explicitly through the \code{formula} argument.
 
If \code{estmethod="ts"}, then two-stage estimation of \eqn{\psi} is performed.
In this case, the model \code{fitX.LZ} is used to construct predictions
\eqn{\hat{X}=\hat{E}(X|L,Z)}. These predictions are subsequently used to re-fit 
the model \code{fitY.LX}, with \eqn{X} replaced with \eqn{\hat{X}}. The obtained 
coefficient(s) for \eqn{\hat{X}} in the re-fitted model is the two-stage estimator of \eqn{\psi}. 

If \code{estmethod="g"}, then G-estimation of \eqn{\psi} is performed. In this case,
the estimator is obtained as the solution to the estimating equation
\deqn{H(\psi)=\sum_{i=1}^n\hat{d}(L_i,Z_i)h_i(\psi)=0.}
The function \eqn{h_i(\psi)} is defined as
\deqn{h_i(\psi)=Y_i-m^T(L_i)\psi X_i}
when \code{link="identity"},
\deqn{h_i(\psi)=Y_i\textrm{exp}\{-m^T(L_i)\psi X_i\}}
when \code{link="log"}, and
\deqn{h_i(\psi)=\textrm{expit}[\textrm{logit}\{\hat{E}(Y|L_i,Z_i,X_i)\}-m^T(L_i)\psi X_i]}
when \code{link="logit"}. In the latter, \eqn{\hat{E}(Y|L_i,Z_i,X_i)} is
an estimate of \eqn{E(Y|L_i,Z_i,X_i)} obtained from the model \code{fitY.LZX}.
The estimated function \eqn{\hat{d}(L,Z)} is chosen so that the true function
has conditional mean 0, given \eqn{L}; \eqn{E\{d(L,Z)|L\}=0}. 
The specific form of \eqn{\hat{d}(L,Z)} is determined by the user-specified models. 
If \code{fitX.LZ} and \code{fitX.L} are specified, then \eqn{\hat{d}(L,Z)=m(L)\{\hat{E}(X|L,Z)-\hat{E}(X|L)\}}, 
where \eqn{\hat{E}(X|L,Z)} and \eqn{\hat{E}(X|L)} are obtained from \code{fitX.LZ}
and \code{fitX.L}, respectively. If these are not specified, then \eqn{\hat{d}(L,Z)=m(L)\{Z-\hat{E}(Z|L)\}}, 
where \eqn{\hat{E}(Z|L)} is obtained from \code{fitZ.L}, which then must be specified. 
}

\value{
\code{ivglm} returns an object of class \code{"ivglm"}, which inherits from 
class \code{"ivmod"}. An object of class \code{"ivglm"} is a list containing 
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
\code{fitY.LX}. If \code{estmethod="g"}, then the first columns correspond to \eqn{\psi}, 
and the remaining columns correspond to the parameters estimated by \code{fitZ.L},
\code{fitX.LZ}, \code{fitX.L} and \code{fitY.LZX}, whichever were used in the
estimation process.    
  }
\item{d.estfun}{
  the jacobian matrix of \code{colMeans(estfunall)}.    
  }
\item{converged}{
  logical. Was a solution found to the estimating equations?   
}
\item{fitY.LX}{
  the re-fitted model \code{fitY.LX} used in the estimation process when \code{estmethod="ts"}.
  This element is NULL when \code{estmethod="g"}. 
}

}

\references{
Bowden J., Vansteelandt S. (2011). Mendelian randomization analysis of
case-control data using structural mean models. \emph{Statistics in Medicine}
\bold{30}(6), 678-694.

Sjolander A., Martinussen T. (2019). Instrumental variable estimation with the R package ivtools.
\emph{Epidemiologic Methods} \bold{8}(1), 1-20.
 
Vansteelandt S., Bowden J., Babanezhad M., Goetghebeur E. (2011). 
On instrumental variables estimation of causal odds ratios.
\emph{Statistical Science} \bold{26}(3), 403-422.

}

\author{
Arvid Sjolander. 
}

\note{
\code{ivglm} allows for weights. However, these are defined implicitly 
through the input models. Thus, when models are used as input to \code{ivglm}, 
these models have to be fitted with the same weights. When \code{estmethod="g"}
the weights are taken from \code{fitX.LZ}, if specified by the user. If \code{fitX.LZ} is not
specified then the weights are taken from \code{fitZ.L}. Hence, if weights are used,
then either \code{fitX.LZ} or \code{fitZ.L} must be specified.
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

#two-stage estimation
fitX.LZ <- glm(formula=X~Z, data=data)
fitY.LX <- glm(formula=Y~X+L+X*L, data=data)
fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=data, 
  ctrl=TRUE) 
summary(fitIV)

#G-estimation with model for Z
fitZ.L <- glm(formula=Z~L, data=data)
fitIV <- ivglm(estmethod="g", X="X", Y="Y", fitZ.L=fitZ.L, data=data, 
  formula=~L, link="identity")
summary(fitIV)

#G-estimation with model for X
fitX.LZ <- glm(formula=X~L+Z, data=data)
fitX.L <- glm(formula=X~L, data=data)
fitIV <- ivglm(estmethod="g", Y="Y", fitX.LZ=fitX.LZ, fitX.L=fitX.L, data=data, 
  formula=~L, link="identity")
summary(fitIV)

##---Example 2: logistic model and no covariates--- 

Z <- rbinom(n, 1, 0.5)
X <- rbinom(n, 1, 0.7*Z+0.2*(1-Z)) 
m0 <- plogis(1+0.8*X-0.39*Z)
Y <- rbinom(n, 1, plogis(psi0*X+log(m0/(1-m0)))) 
data <- data.frame(Z, X, Y)

#two-stage estimation
fitX.LZ <- glm(formula=X~Z, family="binomial", data=data)
fitY.LX <- glm(formula=Y~X, family="binomial", data=data)
fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=data, 
  ctrl=TRUE) 
summary(fitIV)

#G-estimation with model for Z
fitZ.L <- glm(formula=Z~1, data=data)
fitY.LZX <- glm(formula=Y~X+Z+X*Z, family="binomial", data=data)
fitIV <- ivglm(estmethod="g", X="X", fitZ.L=fitZ.L, fitY.LZX=fitY.LZX, 
  data=data, link="logit")
summary(fitIV)


}
