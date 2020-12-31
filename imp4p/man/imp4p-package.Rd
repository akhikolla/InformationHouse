\name{imp4p-package}
\alias{imp4p-package}
\alias{imp4p}
\docType{package}
\title{Introduction to the IMP4P package}
\description{
This package provides functions to analyse missing value mechanisms in the context of bottom-up MS-based quantitative proteomics.

It allows estimating a mixture model of missing completely-at-random (MCAR) values and missing not-at-random (MNAR) values.

It also contains functions allowing the imputation of missing values under hypotheses of MCAR and/or MNAR values.

The main functions of this package are the \code{\link{estim.mix}} (estimation of a model of MCAR and MNAR (left-censored) values),  \code{\link{impute.mi}} (multiple imputation) and \code{\link{impute.mix}} (imputation based on a decision rule). It provides also several imputation algorithms for MS-based data. They can be used to impute matrices containing peptide intensities (as Maxquant outputs for instance).

Missing values has to be indicated with NA and a log-2 transformation of the intensities has to be applied before using these functions.

More explanations and details on the functions of this package are available in Giai Gianetto et al.(2020).
}
%\details{
%The DESCRIPTION file:
%\packageDESCRIPTION{imp4p}
%\packageIndices{imp4p}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
%}
\author{
%\packageAuthor{imp4p}

Maintainer: Quentin Giai Gianetto <quentin2g@yahoo.fr>%\packageMaintainer{imp4p}
}
\references{
Giai Gianetto, Q., Wieczorek S., Couté Y., Burger, T. (2020).  A peptide-level multiple imputation strategy accounting for the different natures of missing values in proteomics data. bioRxiv 2020.05.29.122770; doi: \href{https://doi.org/10.1101/2020.05.29.122770}{10.1101/2020.05.29.122770}
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{}
