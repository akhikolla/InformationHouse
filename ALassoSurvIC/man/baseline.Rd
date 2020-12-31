\name{baseline}
\alias{baseline}
\alias{baseline.default}
\alias{print.baseline}
\title{Obtaining the nonparametric maximum likelihood estimate (NPMLE) for the baseline cumulative hazard function}
\usage{
\method{baseline}{default}(object, ...)
}

\arguments{
\item{...}{for S4 method only.}
\item{object}{the \code{object} must be the object retruned by the \code{alacoxIC} function or the \code{unpencoxIC} function.}
}

\description{
Extracting the NPMLE for the baseline cumulative hazard function from the input object. The input object must be the objects returned by the \code{alacoxIC} function or the \code{unpencoxIC} function. The support set over which the cumulative hazard increases is the same as that of the nonparametric maximum likelihood estimator, characterized by Alioum and Commenges (1996). The full details are available in Li et al. (2019).
}

\details{
The estimator for the baseline cumulative hazard function increases only on some support sets, so called maximal intersections, and the NPMLE is indifferent to how it increases on the support sets. The definition of maximal intersections and other details are available in Alioum and Commenges (1996) and Li et al. (2019).
}

\value{
A list with components:
\item{support}{
The maximal intersections with a finite upper endpoint.
}
\item{lambda}{
The jump sizes over the support set.
}
\item{cum.lambda}{
The NPMLE of the baseline cumulative hazard function.
}
}


\references{Alioum, A. and Commenges, D. (1996). A proportional hazards model for arbitrarily censored and truncated data. Biometrics 52, 512-524.

Li, C., Pak, D., & Todem, D. (2019). Adaptive lasso for the Cox regression with interval censored and possibly left truncated data. Statistical methods in medical research. \href{https://doi.org/10.1177/0962280219856238}{https://doi.org/10.1177/0962280219856238}}

\examples{
library(ALassoSurvIC)

### Display the hazard function for the interval censored data
data(ex_ICLT) # the 'virtual' data having 100 subjects and 6 covariates
lowerIC <- ex_ICLT$lowerIC
upperIC <- ex_ICLT$upperIC
trunc <- ex_ICLT$trunc
X <- ex_ICLT[, -c(1:3)]
result <- unpencoxIC(lowerIC, upperIC, X, trunc)
baseline(result)
}

\seealso{\code{alacoxIC}; \code{unpencoxIC}}
\keyword{methods}
