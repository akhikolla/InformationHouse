\name{HWPerm.mult}
\alias{HWPerm.mult}
\title{
Permutation tests fo variants with multiple alleles
}
\description{
  Function \code{HWPerm.mult} implements permutation tests for
  Hardy-Weinberg equilibrium for autosomal and X-chromosomal variants.
}
\usage{
HWPerm.mult(x, y = NULL, nperm = 17000, eps = 1e-10, verbose = TRUE, ...)
}
\arguments{
  \item{x}{vector or triangular matrix with male genotype counts}
  \item{y}{vector or triangular matrix with female genotype counts}
  \item{nperm}{number of permutations (17.000 by default)}
  \item{eps}{a tolerance for the comparison of floating point numbers}
  \item{verbose}{print output or not}
  \item{\dots}{addtional arguments}
}
\details{
  This function approximates exact test probabilities for joint tests
  for HWE and equality of allele frequencies for variants with multiple
  alleles. For purely bi-allelic variant \code{HWPerm} can be used which
  allows for more statistics than just probabilities.

  If argument \code{y} is not specified, gender is considered
  irrelevant, and \code{x} contains total genotype counts. If \code{x}
  and \code{y} are specified, \code{x} should contain male genotype
  counts and \code{y} female genotype counts. \code{x} and \code{y} can
  be vectors if the variant is bi-allelic, but are assumed lower
  triangular if there are more than two alleles. \code{x} is still a
  vector if there are multiple alleles but the variant is
  X-chromosomal. See the examples given below.
}
\value{
  \item{pofthesample}{probability of the observed sample}
  \item{pseudodist}{probabilities of simulated samples}
  \item{pval}{p-value}
}
\references{
  Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that account for gender.
  \href{https://doi.org/10.1111/1755-0998.12748}{https://doi.org/10.1111/1755-0998.12748}
}
\author{
  Jan Graffelman \email{jan.graffelman@upc.edu}
}
\seealso{
  \code{\link{HWPerm}}
}
\examples{
#
# bi-allelic autosomal
#

x1 <- c(AA=298,AB=489,BB=213)
\dontrun{
out <- HWPerm.mult(x1)
}

#
# bi-allelic X-chromosomal
#

x2.m <- c(A=39, B=21)
x2.f <- toTriangular(c(AA=28, AB=30, BB=8))
\dontrun{
out <- HWPerm.mult(x2.m,x2.f)
}

#
# autosomal k alleles not accounting for gender
#

x3 <- c(AA=12,AB=19,AC=13,BB=7,BC=5,CC=0)
x3 <- toTriangular(x3)
\dontrun{
out <- HWPerm.mult(x3)
}

#
# X-chromosomal k alleles 
#

x4.m <- c(A=15,B=17,C=24) 
x4.f <- toTriangular(c(AA=4,AB=2,AC=13,BB=6,BC=19,CC=4))
\dontrun{
out <- HWPerm.mult(x4.m,x4.f)
}

#
# Autosomal k alleles accounting for gender
#

x5.m <- toTriangular(c(AA=12,AB=19,AC=13,BB=7,BC=5,CC=0))
x5.f <- toTriangular(c(AA=8,AB=12,AC=13,BB=8,BC=7,CC=0))
\dontrun{
out <- HWPerm.mult(x5.m,x5.f)
}

}
\keyword{misc}

