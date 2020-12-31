\name{cgAUC-package}
\alias{cgAUC-package}
\docType{package}
\title{
Calculate AUC when gold standard is continuous with large variables.
}
\description{
In this package, the cgAUC is used to calculate the AUC-type measure raised in Obuchowski(2006) when gold standard is continuous.
}
\details{
\tabular{ll}{
Package: \tab cgAUC\cr
Type: \tab Package\cr
Version: \tab 1.2.1\cr
Date: \tab 2014-08-24\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Yuan-chin I. Chang, Yu-chia Chang, and Ling-wan Chen

Maintainer: Yu-chia Chang <curare7177@gmail.com>
}
\references{
Chang, YCI. Maximizing an ROC type measure via linear combination of markers when the gold reference is continuous. Statistics in Medicine 2012.\cr
Obuchowski NA. An ROC-type measure of diagnostic accuracy when the gold standard is continuous-scale. Statistics in Medicine 2006; 25:481--493.\cr
Obuchowski N. Estimating and comparing diagnostic tests accuracy when the gold standard is not binary. Statistics in Medicine 2005; 20:3261--3278.\cr
Friedman JH, Popescu BE. Gradient directed regularization for linear regression and classification. Technical Report, Department of Statistics, Stanford University, 2004.\cr
}
\keyword{ cgAUC }
\examples{
%% ~~ simple examples of the most important functions ~~
# n = 100; p = 5;
# r.x = matrix(rnorm(n * p), , p) # raw data
# r.z = r.x[ ,1] + rnorm(n) # gold standard
# x = scale(r.x) # standardized of raw data
# z = scale(r.z) # standardized of gold standard
# h = n^(-1 / 2)
# t1 = cgAUC(r.x, r.z, h, delta = 1, auto = FALSE, tau = 1, scale = 1) # the delta be constant
# t1
# t2 = cgAUC(r.x, r.z, h, delta = 1, auto = TRUE, tau = 1, scale = 1) # the delta be variable
# t2
}
