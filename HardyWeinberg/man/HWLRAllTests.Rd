\name{HWLRAllTests}
\alias{HWLRAllTests}
\title{
Perform most relevant likelihood ratio test for Hardy-Weinberg
equilibrium and equality of allele frequencies
}
\description{
Function \code{HWLRAllTests} performs a set of likelihood ratio tests in
relation with Hardy-Weinberg proportions (HWP) and equality of allele
frequencies (EAF) for autosomal bi-allelic genetic variants.
}
\usage{
HWLRAllTests(x, y)
}
\arguments{
  \item{x}{Male genotype counts (AA,AB,BB)}
  \item{y}{Female genotype counts (AA,AB,BB)}
}
\details{
Function \code{HWLRAllTests} calls \code{HWLRtest} and calculates the
p-value of six different tests: 1) joint HWP and EAF (A-F); 2) EAF
irrespective of HWP (C-F); 3) HWP irrespective of EAF (D-F); 4) HWP
versus EIC (given EAF) (A-B); 5) EIC irrespective of EAF (E-F) and 6)
HWP versus EIC. Letters refer to scenarios described by Graffelman &
Weir (2017).
}
\value{A named vector with six p-values}
\references{
 Graffelman, J. and Weir, B.S. (2017) On the testing of Hardy-Weinberg
 proportions and equality of allele frequencies in males and females at
 bi-allelic genetic
 markers. \href{http://dx.doi.org/10.1002/gepi.22079}{Genetic Epidemiology} pp. 1--15.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu}}
\seealso{
  \code{\link{HWLRtest}}
  }
\examples{
males <- c(AA=11,AB=32,BB=13) 
females <- c(AA=14,AB=23,BB=11)
pvalues <- HWLRAllTests(males,females)
print(pvalues)
}
\keyword{misc}
