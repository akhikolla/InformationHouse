\name{linearQ-package}
\alias{linearQ-package}
\alias{linearQ}
\docType{package}
\title{
   Quantile Simulation in Multiscale Chang-Point Segmentation.
}
\description{
  A linear algorithm for simulating quantiles of multiscale statistics under null hypothesis in multiscale change-point segmentation.
}
\details{
  \tabular{ll}{
Package: \tab linearQ \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-02-27\cr
License: \tab GPL (>= 2)\cr
}
Index:
\preformatted{
fastQuantile      Quantile simulation
}
}
\author{
\packageAuthor{linearQ}

Maintainer: \packageMaintainer{linearQ}
}
\references{
  Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference. 
\emph{J. R. Statist. Soc. B, with discussion and rejoinder by the authors}, 76:495--580.
  Li, H., Munk, A., and Sieling, H. (2015). FDR-control in multiscale change-point segmentation. arXiv:1412.5844.
}
\keyword{ package }

\examples{
  # simulate quantiles for multiscale statistics from Normal regression model
    seed = 123
    q  <- fastQuantile(0.9, 500, 100, seed = seed, type = 1)
}
