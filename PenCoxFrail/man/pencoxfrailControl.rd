\name{pencoxfrailControl}
\alias{pencoxfrailControl}
\concept{pencoxfrailControl}
\title{Control Values for \code{pencoxfrail} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{pencoxfrail} function.
}

\usage{

pencoxfrailControl(start = NULL, q_start = NULL, conv.eps = 1e-4, 
                          standardize = FALSE, center = FALSE,
                          smooth=list(nbasis = 6, penal = 0.1), 
                          ridge.pen = 1e-4, print.iter = FALSE, 
                          max.iter = 100, c.app = 1e-6, zeta = 0.5, 
                          exact = 1e-2, xr = NULL, ...)
} 
    
\arguments{
  \item{start}{a vector of suitable length containing starting values for the spline-coefficients of the baseline hazard and the time-varying effects, followed by the fixed and random effects. The correct ordering is important. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal, depending on the dimension of the random effects.}
  \item{conv.eps}{controls the speed of convergence. Default is 1e-4.}
  \item{center}{logical. If true, the covariates corresponding to the time-varying effects will be
    centered. Default is FALSE (and centering is only recommended if really necessary; it can also have a strong effect on the baseline hazard, in particular, if a strong penalty is selected).}
  \item{standardize}{logical. If true, the the covariates corresponding to the time-varying effects will be
    scaled to a variance equal to one (*after* possible centering). Default is FALSE.}
      \item{smooth}{a list specifying the number of basis functions \code{nbasis} (used for the baseline hazard and all time-varying effects) and the smoothness penalty parameter \code{penal}, which is only applied to the baseline hazard. All time-varying effects are penalized by the specific double-penalty \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)} (see \code{\link{pencoxfrail}}), which is based on the overall penalty parameter \eqn{\xi} (specified in the main function \code{\link{pencoxfrail}}) and on the weighting between the two penalty parts \eqn{\zeta}. The degree of the B-splines is fixed to be three (i.e. cubic splines).}
  \item{ridge.pen}{On all time-varying effects (except for the baseline hazard) a slight ridge penalty is applied on the second order differences of the corresponding spline coefficients to stabilize estimation. Default is 1e-4.}
    \item{print.iter}{logical. Should the number of iterations be printed? Default is FALSE.}
  \item{max.iter}{the number of iterations for the final Fisher scoring re-estimation procedure. Default is 200.}
    \item{c.app}{The parameter controlling the exactness of the quadratic approximations of the penalties. Default is 1e-6.}
      \item{zeta}{The parameter controlling the weighting between the two penalty parts in the specific double-penalty \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)} (see \code{\link{pencoxfrail}}). Default is 0.5.}
\item{exact}{controls the exactness of the (Riemann) integral approximations. Default is 1e-2.}
\item{xr}{maximal time point that is regarded. Default is NULL and the maximal event or censoring time point in the data is used.}
\item{...}{Futher arguments to be passed.}
}

\value{
  a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{groll@math.lmu.de}
}

\seealso{
  \code{\link{pencoxfrail}}
}

\examples{
# Use different weighting of the two penalty parts
# and lighten the convergence criterion 
pencoxfrailControl(zeta=0.3, conv.eps=1e-3)
}
