\name{augSIMEX-package}
\alias{augSIMEX-package}
\docType{package}
\title{
  R package for Analysis of Data with Mixed Measurement Error and Misclassification in Covariates
}
\description{
Implementation of the augmented SIMEX algorithm for data with mixed measurement error and misclassification in covariates. The package is allowing for both parametric SIMEX and empirical SIMEX.
}
\details{
  The \code{SIMEX} package implements the method developed by Yi et al. (2015) to adjust for both measurement error and misclassification in the generalized linear models. The main function is \code{\link[augSIMEX]{augSIMEX}} which returns an "augSIMEX" object. The user can summarize the returned augSIMEX object in a similar format as in \code{\link[stats]{glm}}. The extrapolation can also be visualized through \code{plot} function. Other implementation tools are also provided in the package.
}
\author{
Qihuang Zhang and Grace Y. Yi.

Maintainer: Qihuang Zhang <qihuang.zhang@uwaterloo.ca>
}
\references{
  Yi G Y, Ma Y, Spiegelman D, et al. Functional and structural methods with mixed measurement error and misclassification in covariates[J]. Journal of the American Statistical Association, 2015, 110(510): 681-696.
}
\keyword{ package }
\seealso{
  \code{\link[augSIMEX]{augSIMEX}}
}

