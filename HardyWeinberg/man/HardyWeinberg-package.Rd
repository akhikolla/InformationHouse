\name{HardyWeinberg-package}
\alias{HardyWeinberg-package}
\alias{HardyWeinberg}
\docType{package}
\title{Statistical Tests and Graphics for Hardy-Weinberg Equilibrium}
\description{
   Contains tools for exploring Hardy-Weinberg equilibrium (Hardy, 1908;  Weinberg, 1908) for bi and multi-allelic genetic marker data. All classical tests (chi-square, exact, likelihood-ratio and permutation tests) with bi-allelic variants are included in the package, as well as functions for power computation and for the simulation of marker data under equilibrium and disequilibrium. Routines for dealing with markers on the X-chromosome are included (Graffelman & Weir, 2016) including Bayesian procedures. Some exact and permutation procedures also work with multi-allelic variants. Special test procedures that jointly address Hardy-Weinberg equilibrium and equality of allele frequencies in both sexes are supplied, for the bi and multi-allelic case. Functions for testing equilibrium in the presence of missing data by using multiple imputation are also provided. Implements several graphics for exploring the equilibrium status of a large set of bi-allelic markers: ternary plots with acceptance regions, log-ratio plots and Q-Q plots. 
}
\details{
\tabular{ll}{
Package: \tab HardyWeinberg\cr
Type: \tab Package\cr
Version: \tab 1.7.0\cr
Date: \tab 2020-12-13\cr
License: \tab GPL Version 2 or later.\cr
}

With function \code{HWTernaryPlot} one can create ternary plots with
acceptance regions for HWE. Several routines implement statistical tests
for HWE such as \code{HWChisq}, \code{HWExact}, \code{HWLratio} and
\code{HWPerm}. Bayesian procedures are available using
\code{HWPosterior}. Aikaike information criterion for various scenarios
can be calculated with \code{HWAIC}. Power computations are possible
with \code{HWPower}. 

}
\author{
Jan Graffelman

Maintainer: Jan Graffelman <jan.graffelman@upc.edu>

}
\references{
Weir, B.S. (1996) \emph{Genetic Data Analysis II}. Sinauer Associates, Massachusetts.

Graffelman, J. and Morales, J. (2008) Graphical Tests for Hardy-Weinberg Equilibrium 
Based on the Ternary Plot. \emph{Human Heredity} 65(2):77-84. 

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The {HardyWeinberg} Package.
\emph{Journal of Statistical Software} 64(3):
1-23. \url{https://www.jstatsoft.org/v64/i03/}.

Graffelman, J. and Weir, B.S. (2016) Testing for Hardy-Weinberg
equilibrium at bi-allelic genetic markers on
the X chromosome. \emph{Heredity} 116(6):
558-568. \url{https://doi.org/10.1038/hdy.2016.20} Nature Publishing Group.

Hardy, G.H. (1908) Mendelian Proportions in a Mixed Population.
\emph{Science} 28(706): 49-50. 

Weinberg, W. (1908) On the Demonstration of Heredity in Man. In: Papers on Human Genetics.
S. H. Boyer (Ed.) Englewood Cliffs, NJ:Prentice Hall, translated, 1963.
}
\keyword{ package }
\examples{

library(HardyWeinberg)

# draw random SNPs from a population that is in HWE

set.seed(123)

m <- 100 # number of markers
n <- 100 # sample size

X <- HWData(n,m)
out <- HWTernaryPlot(X,100,region=1,vertex.cex=2,signifcolour=TRUE)

}
