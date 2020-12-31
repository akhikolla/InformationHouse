%% Define Macros to copy information from DESCRIPTION file
% To get the package title at build time from the DESCRIPTION file
\newcommand{\packageTitle}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_title("#1")}}
% To get the package description at build time from the DESCRIPTION file
\newcommand{\packageDescription}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_description("#1")}}
% To get the package author at build time from the DESCRIPTION file
\newcommand{\packageAuthor}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_author("#1")}}
% To get the package maintainer at build time from the DESCRIPTION file
\newcommand{\packageMaintainer}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer("#1")}}
% To get a formatted copy of the whole DESCRIPTION file
\newcommand{\packageDESCRIPTION}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_DESCRIPTION("#1")}}
% To include various indices about an installed package
\newcommand{\packageIndices}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_indices("#1")}}

\name{dHSIC-package}
\alias{dHSIC-package}
\alias{dHSIC}
\docType{package}
\title{\packageTitle{dHSIC}}
\description{\packageDescription{dHSIC}}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{dHSIC}
\packageIndices{dHSIC}
}
\author{
\packageAuthor{dHSIC}

Maintainer: \packageMaintainer{dHSIC}
}
\references{
  Gretton, A., K. Fukumizu, C. H. Teo, L. Song, B. Sch{\"o}lkopf and
  A. J. Smola (2007). A kernel statistical test of independence. In
  Advances in Neural Information Processing Systems (pp. 585-592).
  
  Pfister, N., P. B{\"u}hlmann, B. Sch{\"o}lkopf and J. Peters (2017).
  Kernel-based Tests for Joint Independence. To appear in the
  Journal of the Royal Statistical Society, Series B.
}

\keyword{ package }
