\name{int.approx}
\alias{int.approx}
\concept{int.approx}
\title{Approximation of a Cox likelihood intergral}
\description{The function approximates the integral \eqn{\int_0^t exp(u B(s) \alpha) ds)} which appears in the (full) Cox likelihood
if the covariate \eqn{u} has a time-varying effect \eqn{\beta(t)}, which is expanded in B-splines, i.e. \eqn{\beta(t) = B(t) \alpha}.
}

\usage{
int.approx(z,time.grid,B,nbasis,alpha)
} 
    
\arguments{
  \item{z}{a vector which contains at the first component a time point up to which it should be integrated and the covariates \eqn{u} in the remaining components.}
  \item{time.grid}{an equally-spaced time grid on which the B-spline design matrix \eqn{B} has been generated.
  The maximal value of the time grid should usually be the maximal upper integral border that is of interest.}
  \item{B}{a B-spline design matrix, which has been created with the function \code{\link{bs.design}} on the full time grid \code{time.grid}.}
  \item{nbasis}{number of basis functions used when the B-spline design matrix \eqn{B} has been generated.}
  \item{alpha}{vector of B-spline coefficients.}
}

\value{
  The B-spline design matrix is returned.
}

\author{
Andreas Groll \email{groll@math.lmu.de}
}

\seealso{
  \code{\link{pencoxfrail}}
}

\examples{
## generate time grid and corresponding B-spline design matrix
time.grid <- seq(0,200,by=1)
B <- bs.design(x=time.grid, xl=min(time.grid), xr=max(time.grid), spline.degree=3, nbasis=5)

## specify spline coefficients and covariate vector (with upper integral bound as first component)
alpha <- c(0.1,0.2,0.05,0.1,0.15)
z <- c(time=100,age=25)

## calculate intergal from 0 to 100
int.approx(z=z,time.grid=time.grid,B=B,nbasis=5,alpha=alpha)
}
