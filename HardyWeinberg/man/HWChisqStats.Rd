\name{HWChisqStats}
\alias{HWChisqStats}
\title{
Fast computation of chi-square statistics for Hardy-Weinberg equilibrium
}
\description{
\code{HWChisqStats} is a function for the fast computation of chi-square
statistics (or the corresponding p-values) for a large set of bi-allelic markers (typically SNPs). 
}
\usage{
HWChisqStats(X, x.linked = FALSE, pvalues = FALSE)
}
\arguments{
  \item{X}{A matrix with genotype counts, one row per marker. \code{X}
    should have 5 columns for an X-chromosomal data set and 3 columns
    for an autosomal data set.}
  \item{x.linked}{Logical indicating whether the markers are autosomal (\code{x.linked=FALSE})
    or X-chromosomal (\code{x.linked=TRUE}).}
  \item{pvalues}{Logical indicated whether chi-square statistics should be returned (\code{pvalues=FALSE}) or whether p-values should be returned (\code{pvalues=TRUE}).}
}
\details{
Matrix \code{X} should strictly comply with the following format. For
an autosomal dataset is should contain the 3 genotype counts in order
(AA,AB,BB). For an X-chromosomal dataset it should contain the 5
genotype counts in order (A,B,AA,AB,BB) where A and B are the male
counts and AA, AB and BB the female counts.

This function was written for speed improvement, and should be much
faster than looping over the rows of \code{X} with \code{HWChisq}. There
is no error checking on the supplied data matrix.
}
\value{A vector of chi-square statistics}
\references{
  Graffelman, J. and Weir, B.S. (2016) Testing for Hardy-Weinberg
  equilibrium at 2 bi-allelic genetic markers on the X chromosome.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
#
# Autosomal example
#
set.seed(123)
X <- HWData(1000,100)
monom <- (X[,2]==0 & X[,1]==0) | (X[,2]==0 & X[,3]==0)
X <- X[!monom,] # exclude monomorphics
Chisq.stats <- HWChisqStats(X,x.linked=FALSE,pvalues=FALSE)
Chisq.pvals <- HWChisqStats(X,x.linked=FALSE,pvalues=TRUE)
#
# X-chromosomal example
#
X <- HWData(1000,100,n.males=50,nA=75,x.linked=TRUE)
Chisq.stats <- HWChisqStats(X,x.linked=TRUE,pvalues=FALSE)
Chisq.pvals <- HWChisqStats(X,x.linked=TRUE,pvalues=TRUE)
}
\keyword{htest}
