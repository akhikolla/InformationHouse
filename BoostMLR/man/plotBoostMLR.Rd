\name{plotBoostMLR}
\alias{plotBoostMLR}
\title{
Plotting results across across the boosting iterations.
}
\description{
Plotting training and test error, and estimate of variance/correlation parameters
across the boosting iterations.
}
\usage{
plotBoostMLR(Result,
             xlab = "",
             ylab = "",
             legend_fraction_x = 0.10,
             legend_fraction_y = 0,
             ...)
}
\arguments{
  \item{Result}{
Result in the matrix form either training or test error, or estimate of variance/correlation
parameters across the boosting iterations.
}
\item{xlab}{
Label for the x-axis.
}
  \item{ylab}{
Label for the y-axis.
}
  \item{legend_fraction_x}{
Value use to expland the x-axis.  
}
  \item{legend_fraction_y}{
Value use to expland the y-axis.}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
Plotting training and test error, and estimate of variance/correlation parameters
across the boosting iterations.
}
\author{
  Amol Pande and Hemant Ishwaran
}
\examples{
\donttest{
##-----------------------------------------------------------------
## Multivariate Longitudinal Response
##-----------------------------------------------------------------

# Simulate data involves 3 response and 4 covariates

dta <- simLong(n = 100, N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                            y = dta$y, M = 100, VarFlag = FALSE)

# Plot training error
plotBoostMLR(boost.grow$Error_Rate,xlab = "m",ylab = "Training Error")
}
}
\keyword{plot}
