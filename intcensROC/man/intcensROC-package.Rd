\name{intcensROC-package}
\alias{intcensROC-package}
\alias{intcensROC-package}
\docType{package}
\title{
Fast Spline Function Based Constrained Maximum Likelihood Estimator for 
ROC Estimation of Interval Censored Survival Data
}
\description{
The kernel of this \code{Rcpp} based package is a efficient implementation 
of the generalized gradient projection method for spline function based 
constrained maximum likelihood estimator (Wu, Yuan; Zhang, Ying. Partially 
monotone tensor spline estimation of the joint distribution function 
with bivariate current status data. Ann. Statist. 40, 2012, 1609-1636 
<doi:10.1214/12-AOS1016>). The key function is to compute density 
function of the joint distribution of event time and the marker. 
The core function returns the receiver operating characteristic (ROC) curve for the interval 
censored survival data as well as area under the curve (AUC). 
}
\details{
\tabular{ll}{
Package: \tab intcensROC \cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2020-11-11\cr
License: \tab GPL-3\cr
}
Please refer to the individual 
function documentation or the included vignette for more information. 
The package vignette serves as a tutorial for using this package. 
The technical details are provided in the reference cited below.
}

\author{
Jiaxing Lin, Yuan Wu, Xiaofei Wang, Kouros Owzar.
Maintainer: Jiaxing Lin<jiaxing.lin@duke.edu>
}
\references{
Wu, Yuan; Zhang, Ying. Partially monotone tensor spline estimation of the 
joint distribution function with bivariate current status data. 
Ann. Statist. 40, 2012, 1609-1636 <doi:10.1214/12-AOS1016>
}

