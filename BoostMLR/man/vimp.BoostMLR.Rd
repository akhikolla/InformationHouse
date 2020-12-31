\name{vimp.BoostMLR}
\alias{vimp.BoostMLR}
\title{Variable Importance}
\description{
  Calculate standardized variable importance (VIMP) for each covariate or a joint VIMP
  of multiple covariates.
}
\usage{
\method{vimp}{BoostMLR}(Object,
     xvar.names = NULL,
     joint = FALSE)}
\arguments{
  \item{Object}{A boosting object of class \code{(BoostMLR, predict)}.}
  \item{xvar.names}{Names of the x-variables for which VIMP is requested. If NULL, VIMP is calcuated for all the covariates.}
  \item{joint}{Whether to estimate VIMP for each covariate from \code{xvar.names} or a joint VIMP for multiple covariates?}
}
\details{
  Standardized variable importance (VIMP) is calcuated for each covariate or a joint VIMP is calculated for all the covariates specified in \code{xvar.names}.
}
\value{
{If \code{joint} = FALSE, a standardized VIMP for each covariate is obtained otherwisea joint VIMP for all the covariates is obtained. 
The result consists of a list of
length equal to the number of multivariate response.
Each element from the list represents a matrix with number of rows equal to the number of covariates (in case of joint VIMP, the matrix will have a single row) and the number of columns equal to the number of overlapping time intervals + 1 where the first column contains covariate main effects and all other columns contain covariate-time interaction effects.}
}
\author{
  Amol Pande and Hemant Ishwaran
}
\references{
  Pande A., Ishwaran H., Blackstone E.H. (2020).  Boosting for multivariate
  longitudinal response.

  Friedman J.H. Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232, 2001.
}
\examples{
\donttest{
##-----------------------------------------------------------------
## Calculate individual and joint VIMP
##-----------------------------------------------------------------

# Simulate data involves 3 response and 4 covariates

dta <- simLong(n = 100, ntest = 100 ,N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")
dtaL <- dta$dtaL
trn <- dta$trn
# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dtaL$features[trn,], tm = dtaL$time[trn], 
                      id = dtaL$id[trn], y = dtaL$y[trn,], M = 100, VarFlag = FALSE)

boost.pred <- predictBoostMLR(Object = boost.grow, x = dtaL$features[-trn,], 
                               tm = dtaL$time[-trn], id = dtaL$id[-trn], 
                               y = dtaL$y[-trn,], importance = FALSE)
# Individual VIMP                               
Ind_vimp <- vimp.BoostMLR(boost.pred)

# Joint VIMP
Joint_vimp <- vimp.BoostMLR(boost.pred,joint = TRUE)

}
}
\keyword{plot}
