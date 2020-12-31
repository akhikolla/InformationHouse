\name{factorcpt-package}
\alias{factorcpt-package}
\alias{factorcpt}
\docType{package}

\title{
Simultaneous multiple change-point and factor analysis for high-dimensional time series
}

\description{
The package implements a two-stage methodology for consistent multiple change-point detection under factor modelling. It performs multiple change-point analysis on the common and idiosyncratic components separately, and thus automatically identifies their origins. The package also implements the Double CUSUM Binary Segmentation algorithm, which is proposed for multiple change-point detection in high-dimensional panel data with breaks in its mean.
}

\details{
\tabular{ll}{
Package: \tab factorcpt\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2016-12-22\cr
License: \tab GPL (>= 2)\cr
}

The main routines of the package are \code{\link{factor.seg.alg}} and \code{\link{func_dc}}.

}
\author{
Haeran Cho, Matteo Barigozzi, Piotr Fryzlewicz

Maintainer: Haeran Cho <haeran.cho@bristol.ac.uk>
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint. \cr

H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. Electronic Journal of Statistics. 10: 2000-2038.
}

\keyword{change-point, factor model, segmentation}
