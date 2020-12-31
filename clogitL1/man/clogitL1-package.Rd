\name{clogitL1-package}
\alias{clogitL1-package}
\docType{package}
\title{
Penalised conditional logistic regression.
}
\description{
Tools for the fitting and cross validation of exact conditional logistic regression models with lasso and elastic net penalties. Uses cyclic coordinate descent and warm starts to compute the entire path efficiently.
}
\details{
\tabular{ll}{
Package: \tab clogitL1\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2013-05-06\cr
License: \tab GPL-2\cr
}
Very simple to use. The main fitting function \code{clogitL1} accepts x, y data and a strata vector indicating stratum membership. It fits the exact conditional logistic regression model at a grid of regularisation parameters. 
Only 7 functions:
\itemize{
\item \code{clogitL1}
\item \code{cv.clogitL1}
\item \code{plot.clogitL1}
\item \code{plot.cv.clogitL1}
\item \code{print.clogitL1}
\item \code{summary.clogitL1}
\item \code{summary.cv.clogitL1}
}
}

\references{
\url{http://www.jstatsoft.org/v58/i12/}
}

\author{
Stephen Reid and Rob Tibshirani

Maintainer: Stephen Reid <sreid@stanford.edu>
}
