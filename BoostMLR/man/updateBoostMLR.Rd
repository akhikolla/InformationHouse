\name{updateBoostMLR}
\alias{updateBoostMLR}
\title{
Update boosting object with an additional boosting iterations
}
\description{
Function allows to update boosting object with an additional boosting iterations.
}
\usage{
updateBoostMLR(Object,
               M_Add,
               Verbose = TRUE,
               ...)
}
\arguments{
\item{Object}{Boosting object. This object is previously obtained using \command{BoostMLR} function or using \command{update} function.}
\item{M_Add}{Number of additional boosting iterations.}
\item{Verbose}{Print the current stage of boosting iteration?}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
In boosting, \code{Mopt}, the number of boosting iterations
required to achive optimal result, is unknown. 
Typically, \code{Mopt} is estimated by specifying a large value
of \code{M} and then search for an optimal value that is less 
than \code{M} using the test data. 
Function \code{update} allows user to start with a small value of \code{M}, and keep increamenting boosting iterations, each time running through the test data, until an optimal boosting iteration is found. This can significantly reduce unnecessary computations, particularly when \code{Mopt} << \code{M}.
The procedure can be replicated multiple times using the boosting object (see example below).
Results from \command{update} can be treated the same way we treat results
from \command{BoostMLR}.
}
\author{
Amol Pande and Hemant Ishwaran
}

\seealso{
  \command{\link{BoostMLR}},
  \command{\link{predictBoostMLR}},
  \command{\link{simLong}}
}
\examples{
\donttest{
##-----------------------------------------------------------------
## Univariate Longitudinal Response
##-----------------------------------------------------------------

# Simulate data involves 1 response and 4 covariates

dta <- simLong(n = 100, N = 5, rho =.80, model = 2, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                          y = dta$y, M = 100, VarFlag = FALSE)
                          
# Update boosting object for the additional 100 iteration
boost.grow <- updateBoostMLR(Object = boost.grow, M_Add = 100,Verbose = TRUE)

# Update boosting object for the additional 50 iteration
boost.grow <- updateBoostMLR(Object = boost.grow, M_Add = 50,Verbose = TRUE)
}
}
\keyword{boosting}
