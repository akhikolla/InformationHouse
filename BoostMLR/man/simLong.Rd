\name{simLong}
\alias{simLong}
\title{Simulate longitudinal data}
\description{
  Simulates longitudinal data from multivariate and univariate longitudinal 
  response model.}
\usage{simLong(n = 100,
        ntest = 0,
        N = 5,
        rho = 0.8,
        model = c(1, 2),
        phi = 1,
        q_x = 0,
        q_y = 0,
        type = c("corCompSym", "corAR1", "corSymm", "iid"))}
\arguments{
  \item{n}{Requested training sample size.}
  \item{ntest}{Requested test sample size.}
  \item{N}{Parameter controlling number of time points per subject.}
  \item{rho}{Correlation parameter.}
  \item{model}{Requested simulation model.}
  \item{phi}{Variance of measurement error.}
  \item{q_x}{Number of noise covariates.}
  \item{q_y}{Number of noise responses.}
  \item{type}{Type of correlation matrix.}
}
\details{
  Simulates longitudinal data from multivariate and univariate longitudinal 
  response model. We consider following 2 models:
  \enumerate{
    \item \emph{\code{model=1}:} Simpler linear model consist of three
    longitudinal responses, \code{y1}, \code{y2}, and \code{y3} and
    four covariates \code{x1}, \code{x2}, \code{x3}, and \code{x4}. 
    Response \code{y1} is associated with \code{x1} and \code{x4}.
    Response \code{y2} is associated with \code{x2} and \code{x4}.
    Response \code{y3} is associated with \code{x3} and \code{x4}.
    
    \item \emph{\code{model=2}:} Relatively complex model consist of
    single longitudinal response and four covariates. Model includes
    non-linear relationship between response and covariates and
    covariate-time interaction.
  }
}
\value{
  An invisible list with the following components:
  \item{dtaL}{List containing the simulated data in the following order:
    \code{features}, \code{time}, \code{id} and \code{y}.}
  \item{dta}{Simulated data given as a data frame.}
  \item{trn}{Index of \code{id} values identifying the training data.}
}
\author{
  Amol Pande and Hemant Ishwaran
}
\references{
  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
  longitudinal data, \emph{Machine Learning}, 106(2): 277--305. 
}
\keyword{simulation}
\keyword{variable selection}
