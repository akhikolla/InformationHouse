\name{predictBoostMLR}
\alias{predictBoostMLR}
\title{
Prediction for the multivariate longitudinal response
}
\description{
Function returns predicted values for the response. Also, if the response
is provided, function returns test set performance, optimal boosting
iteration, and variable importance (VIMP).
}
\usage{
predictBoostMLR(Object,
                x,
                tm,
                id,
                y,
                M,
                importance = FALSE,
                eps = 1e-5,
                ...)
}
\arguments{
\item{Object}{A boosting object obtained using the function \command{BoostMLR}
on the training data.}
\item{x}{Data frame (or matrix) containing the test set x-values (covariates).   
Covariates can be time-varying or time-invariant.
If \code{x} is unspecified while growing the \code{Object}, it should be
unspecified here as well.}
\item{tm}{Vector of test set time values.
If \code{tm} is unspecified while growing the \code{Object}, it should be
unspecified here as well.}
\item{id}{Vector of test set subject identifier.
If \code{id} is unspecified while growing the \code{Object}, it should be
unspecified here as well.}
\item{y}{Data frame (or matrix) containing the test set y-values 
(response) in case of multivariate response or a 
vector of y-values in case of univariate response.
If \code{y} is unspecified then predicted values corresponding to \code{x}
and \code{tm} can be obtained but no performance measure such as
test set error and VIMP.}
\item{M}{Number of boosting iterations. Value should be less than or equal
to the value specified in the \code{Object}. If unspecified, value from 
the \code{Object} will be used.}
\item{importance}{Whether to calculate standardized variable importance (VIMP) for each covariate?}
\item{eps}{Tolerance value used for determining the optimal \code{M}.}
\item{...}{Further arguments passed to or from other methods.}
}

\details{The predicted response and performance values are obtained for
  the test data using the \code{Object} grown using function \command{BoostMLR} on
  the training data.
}

\value{
\item{Data}{A list with elements \code{x}, \code{tm}, \code{id} 
and \code{y}. Additionally, the list include mean and 
standard deviation of \code{x} and \code{y}.}
\item{x_Names}{Variable names of \code{x}.}
\item{y_Names}{Variable names of \code{y}.}
\item{mu}{Estimate of conditional expectation of \code{y} corresponding
to the last boosting iteration.}
\item{mu_Mopt}{Estimate of conditional expectation of \code{y} corresponding
to the optimal boosting iteration.}
\item{Error_Rate}{Test set error rate for each multivariate response across the boosting iterations.}
\item{Mopt}{The optimal number of boosting iteration.}
\item{nu}{Regularization parameter.}
\item{rmse}{Test set standardized root mean square error (sRMSE) at the \code{Mopt}.}
\item{vimp}{Standardized VIMP for each covariate. This consist of a list of
length equal to the number of multivariate response.
Each element from the list represents a matrix with number of rows equal to the number of covariates and the number of columns equal to the number of overlapping time intervals + 1 where the first column contains covariate main effects and all other columns contain covariate-time interaction effects.}
\item{Pred_Object}{Useful for internal calculation.}
}

\references{
  Pande A., Ishwaran H., Blackstone E.H. (2020).  Boosting for multivariate
  longitudinal response.

  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
  longitudinal data, \emph{Machine Learning}, 106(2): 277--305.
  
  Pande A. (2017).  \emph{Boosting for longitudinal data}.  Ph.D. Dissertation,
  Miller School of Medicine, University of Miami.
}
\author{
Amol Pande and Hemant Ishwaran
}

\seealso{
  \command{\link{BoostMLR}},
  \command{\link{updateBoostMLR}},
  \command{\link{simLong}}
}
\examples{
\donttest{
##-----------------------------------------------------------------
## Multivariate Longitudinal Response
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
                               y = dtaL$y[-trn,], importance = TRUE)
# Plot test set error
plotBoostMLR(boost.pred$Error_Rate,xlab = "m",ylab = "Test Set Error",
                                              legend_fraction_x = 0.2)
}
}
\keyword{boosting}

