\name{ALassoSurvIC-package}
\alias{ALassoSurvIC-package}
\docType{package}
\title{
	Variable selection procedure with the adaptive lasso for interval censored and possibly left truncated data
}
\description{
This package provides penalized variable selection tools for the Cox proportional hazards model with interval censored and possibly left truncated data. The main function \code{alacoxIC} performs the variable selection via a penalized nonparametric maximum likelihood estimation (PNPMLE) with an adaptive lasso penalty. The function also finds the optimal thresholding parameter automatically by minimizing the Bayesian information criterion (BIC). The unpenalized nonparametric maximum likelihood estimation for interval censored and possibly left truncated data is also available with the \code{unpencoxIC} function. The asymptotic validity of the methodology is established in Li et al. (2019).
}

\details{
\tabular{ll}{
Package: \tab ALassoSurvIC\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2019-8-28\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Chenxi Li, Daewoo Pak and David Todem
}
\references{Li, C., Pak, D., & Todem, D. (2019). Adaptive lasso for the Cox regression with interval censored and possibly left truncated data. Statistical methods in medical research. \href{https://doi.org/10.1177/0962280219856238}{https://doi.org/10.1177/0962280219856238}
}
\seealso{
\code{\link{alacoxIC}}; \code{\link{unpencoxIC}}
}
