\name{lclGWAS-package}
\alias{lclGWAS-package}
\alias{lclGWAS}
\docType{package}
\title{
	Efficient Estimation of Multivariate Frailty Model Using Exact Likelihood Function for Grouped Survival Data
}
\description{
The core of this \code{Rcpp} based package is several functions to estimate the baseline hazard, frailty variance, and fixed effect parameter for a discrete-time shared frailty model with random effects. The functions are designed to analyze grouped time-to-event data accounting for family structure of related individuals (i.e., trios). The core functions include two processes: (1) evaluate the multivariable integration to compute the exact proportional hazards model based likelihood and (2) estimate the desired parameters using maximum likelihood estimation. The integration is evaluated by the \code{Cuhre} algorithm from the \code{Cuba} library (Hahn, T., Cuba-a library for multidimensional numerical integration, Comput. Phys. Commun. 168, 2005, 78-95), and the source files of the \code{Cuhre} function are included in this package. The maximization process is carried out using Brent's algorithm, with the \code{C++} code file from John Burkardt and John Denker (Brent, R.,Algorithms for Minimization without Derivatives, Dover, 2002).
License: GPL (>= 2)}
\details{
\tabular{ll}{
Package: \tab lclGWAS\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2017-02-20\cr
License: \tab GPL-3\cr
}
Please refer to the individual function documentation or the included vignette for more information. The package vignette serves as a tutorial for using this package.
}

\author{
Jiaxing Lin, Alexander Sibley, Tracy Truong, Nancy Cox, Eileen Dolan, Yu Jiang, Janice McCarthy, Andrew Allen, Kouros Owzar, Zhiguo Li
Maintainer: Jiaxing Lin <jiaxing.lin@duke.edu>
}

\references{
	Ripatti, S. and Palmgren, J., Estimation of Multivariate Frailty Models Using Penalized Partial Likelihood. Biometrics, 56, 2000, 1016-1022.\cr
	Hahn, T., Cuba-a library for multidimensional numerical integration, Computer Physics Communications, 168, 2005, 78-95.\cr
	Brent, R., Algorithms for Minimization without Derivatives, Dover, 2002.
}

\keyword{ package }
\seealso{
    \code{\link{Rcpp}}
}

