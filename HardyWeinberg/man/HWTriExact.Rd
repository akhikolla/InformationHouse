\name{HWTriExact}
\alias{HWTriExact}
\title{Exact test for Hardy-Weinberg equilibrium and equality of allele
  frequencies for tri-allelic variants.
}
\description{
Function \code{HWTriExact} does a standard exact test for Hardy-Weinberg
equilibrium of a tri-allelic variant, and also does joint exact tests
for equilibrium and equality of allele frequencies if the genotype
counts are given separately for both sexes
}
\usage{
HWTriExact(x, y = NULL, eps = 1e-10, nperm = 17000, verbose = TRUE)
}
\arguments{
  \item{x}{vector with 6 genotype counts (AA,AB,AC,BB,BC,CC)}
  \item{y}{vector with 6 or 3 genotype counts (AA,AB,AC,BB,BC,CC) or (A,B,C)}
  \item{eps}{a tolerance that can be set for comparing exact probabilities}
  \item{nperm}{number of permutations (only relevant for autosomal
  stratified by gender)}
  \item{verbose}{print output or not}
}
\details{
If only \code{x} is specified, an exact test for an autosomal variant
with three alleles will be performed.

If both \code{x} and \code{y} are supplied as vectors with 6 elements,
a permutation test for HWE and equality of allele frequencies (EAF) for an
autosomal variant is performed, using \code{nperm} permutations. The
distribution of the probabilities is returned in \code{pseudodist}. The
computational cost of a completed enumeration algorithm can be
prohibitive in this case.

If \code{x} is supplied as a length 6 vector, and \code{y} as a length 3
vector, the variant is assumed to be X-chromosomal, \code{x} containing
female genotype counts and \code{y} containing male genotyep counts. In
this case a joint exact test for HWE and EAF for an X-chromosomal
tri-allelic variant is executed.

See the examples in the example section below.
}
\value{
  \item{pval}{The p-value of the sample}
  \item{pseudodist}{Distribution of probabilities obtained by simulation}
  \item{pofthesample}{The probability of the observed sample}
}
\references{
  Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that account for gender.
  \href{https://doi.org/10.1111/1755-0998.12748}{https://doi.org/10.1111/1755-0998.12748}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
\code{\link{HWPerm.mult}}
}
\examples{
#
# Autosomal tri-allelic (not accounting for gender)
#

x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
\dontrun{out <- HWTriExact(x)}

#
# Autosomal tri-allelic accounting for gender
#

males   <- c(A=1,B=21,C=34) 
females <- c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15)
\dontrun{out <- HWTriExact(females,males)}

#
# X-chromosomal tri-allelic accounting for gender
#

males   <- c(A=1,B=21,C=34) 
females <- c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15)
\dontrun{out <- HWTriExact(females,males)}
}
\keyword{misc}

