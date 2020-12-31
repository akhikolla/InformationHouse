\name{BoostMLR-package}
\alias{BoostMLR-package}
\docType{package}
\title{
  Boosting for Multivariate Longitudinal Response
}
\description{
The primary feature of the package is to 
jointly model the multiple longitudinal responses (referred to as multivariate longitudinal response) and multiple covariates and time from a longitudinal study using gradient boosting approach (Pande et al., 2020). 
Covariates can be time-varying or time-invariant.
Special cases include modeling of univariate longitudinal
response from a longitudinal study, and univariate or multivariate 
response from a cross-sectional study.
In all cases, responses are assumed to be continuous.
The estimated coefficient can be a function of time (referred to as time-varying coefficient in case of a longitudinal study) or a function of pre-specified covariate (in case of a longitudinal or a cross-sectional study) or fixed.
}
\details{
This package allows joint modeling of a multivariate
longitudinal response, which is based on marginal model.
Estimation is performed using gradient boosting, a generic form
of boosting (Friedman J. H., 2001). The boosting approach use in this package is closely related to
component-wise L2 boosting with L1 penalization.
Package can handle high dimensionalilty of covariate and response
when some of the covariates and responses are pure noise.

The package is designed to identify covariates that affect responses
differently as different time intervals. This idea is helpful
to dissect an overall effect of covariate into different time
intervals. For example, some covariates affect response at
the beginning of the follow-up whereas others at a later stage.
}

\section{Package Overview}{
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

  \enumerate{
    \item \command{\link{BoostMLR}}

    This is the main entry point to the package. The model is fit
    using the gradient boosting approach for the user specified training data.

    \item \command{\link{updateBoostMLR}} (\command{updateBoostMLR})
   
    This allows to update the model by specifying additional boosting iteration.

    \item \command{\link{predictBoostMLR}} (\command{predictBoostMLR})

    Model performace can be obtained using the test set data. This function
    also estimate variable importance (VIMP).

  }
}

\author{
Amol Pande and Hemant Ishwaran

Maintainer: Amol Pande <amoljpande@gmail.com>
}
\references{
  Pande A., Ishwaran H., Blackstone E.H. (2020).  Boosting for multivariate
  longitudinal response.

  Friedman J.H. (2001). Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
}
\keyword{package}
\seealso{
  \command{\link{BoostMLR}},
  \command{\link{updateBoostMLR}},
  \command{\link{predictBoostMLR}},
  \command{\link{simLong}}
}