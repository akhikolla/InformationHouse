\name{HWExactStats}
\alias{HWExactStats}
\title{
Computation of Exact p-values for Hardy-Weinberg equilibrium for sets of SNPs
}
\description{
  \code{HWExactStats} is a function for the computation of Exact p-values
  for a large set of bi-allelic markers (typically SNPs). 
}
\usage{
HWExactStats(X, x.linked = FALSE, plinkcode = TRUE, midp = FALSE,...)
}
\arguments{
  \item{X}{A matrix with genotype counts, one row per marker. \code{X}
    should have 5 columns for an X-chromosomal data set and 3 columns
    for an autosomal data set.}
  \item{x.linked}{Logical indicating whether the markers are autosomal (\code{x.linked=FALSE})
    or X-chromosomal (\code{x.linked=TRUE}).}
  \item{plinkcode}{Logical indicating whether to use faster C++ code
from the PLINK software.}
  \item{midp}{Logical indicating whether to use the mid p-value for the
C++ code or not}
  \item{\dots}{extra arguments that are passed on to \code{HWExact}}
}
\details{
Matrix \code{X} should strictly comply with the following format. For
an autosomal dataset is should contain the 3 genotype counts in order
(AA,AB,BB). For an X-chromosomal dataset it should contain the 5
genotype counts in order (A,B,AA,AB,BB) where A and B are the male
counts and AA, AB and BB the female counts.

Argument \code{plinkcode=TRUE} (the default) will use C++ code for faster calculation (functions \code{SNPHWE2} and
\code{SNPHWEX}) with larger datasets. The C++ code  was generously shared by Christopher
Chang, and the same code is used in the program
\href{http://www.cog-genomics.org/plink/2.0}{PLINK (2.0)}.
}
\value{A vector of p-values}
\references{
  Graffelman, J. and Weir, B.S. (2016) Testing for Hardy-Weinberg
  equilibrium at bi-allelic genetic markers on the X chromosome.
  \href{http://dx.doi.org/10.1038/hdy.2016.20}{Heredity} 116(6) pp. 558--568.
  
  Purcell et al. (2007) PLINK: A Toolset for Whole-Genome Association and
  Population-Based Linkage Analysis. American Journal of Human Genetics
  81(3) pp. 559--575.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} (R code) and
  Christopher Chang \email{chrchang523@gmail.com} (C++ code)}
\seealso{ \code{\link{HWExact}}  }
\examples{
#
# Autosomal example
#
set.seed(123)
X <- HWData(1000,100)
monom <- (X[,2]==0 & X[,1]==0) | (X[,2]==0 & X[,3]==0)
X <- X[!monom,] # exclude monomorphics
Exact.pvalues <- HWExactStats(X,x.linked=FALSE)
#
# X-chromosomal example
#
X <- HWData(1000,100,n.males=50,nA=75,x.linked=TRUE)
Exact.pvalues <- HWExactStats(X,x.linked=TRUE)
}
\keyword{htest}
