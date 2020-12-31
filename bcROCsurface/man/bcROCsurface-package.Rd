\name{bcROCsurface-package}
\alias{bcROCsurface-package}
\alias{bcROCsurface}
\docType{package}
\title{
Bias-Corrected Methods for Estimating the ROC Surface of Continuous Diagnostic Tests
}
\description{
This package provides tools for correcting verification bias in the evaluation of a continuous diagnostic test. More precisely, five bias-corrected methods for ROC surface and VUS inference are provided under MAR assumption, i.e., full imputation (FI), mean score imputation (MSI), inverse probability weighting (IPW), semiparametric efficient (SPE) and K nearest-neighbor (KNN) estimator.
}
\details{
\tabular{ll}{
Package: \tab bcROCsurface\cr
Type: \tab Package\cr
Version: \tab 1.0-4\cr
Date: \tab 2018-11-07\cr
License: \tab GPL 2 | GPL 3\cr
Lazy load: \tab yes\cr
}

Major functions are \code{\link{ROCs}} and \code{\link{vus}}.

}

\author{
Duc-Khanh To, with contributions from Monica Chiogna and Gianfranco Adimari

Maintainer: Duc-Khanh To <toduc@stat.unipd.it>
}

\references{
To Duc, K., Chiogna, M. and Adimari, G. (2018)
Nonparametric estimation of ROC surfaces in presence of verification bias.
\emph{REVSTAT Statistical Journal}. Accepted.

To Duc, K., Chiogna, M. and Adimari, G. (2016)
Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
\emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
}

\keyword{package}

