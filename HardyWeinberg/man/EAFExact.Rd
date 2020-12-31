\name{EAFExact}
\alias{EAFExact}
\title{
Exact test for equality of allele frequencies in males and females
}
\description{
\code{EAFExact} uses a Fisher Exact test to compare allele frequencies
in males and females for variants with k alleles (k >= 2).
}
\usage{
EAFExact(m, f, verbose = TRUE, ...)
}
\arguments{
  \item{m}{vector or triangular matrix with male genotype counts}
  \item{f}{vector or triangular matrix with female genotype counts}
  \item{verbose}{print output (TRUE) or not (FALSE)}
  \item{\dots}{additional arguments for \code{fisher.test}}
}
\details{
  For bi-allelic autosomal variants the genotype counts can be supplied as vectors
  ((AA,AB,BB) for males, and (AA,AB,BB) for females). For X-chromosomal
  bi-allelic variants the genotype counts can also supplied as vectors
  ((A,B) for males, and (AA,AB,BB) for females). For multi-allelic autosomal
  variants male and genotype counts can be supplied as vectors
  (AA,AB,AC,BB,BC,CC,...) or as a triangular
  matrix, where matrix rows and colums are labelled with the allele name
  (A,B,C,...). For multi-allelic X-chromosomal variants, male genotype
  counts must be supplied as a vector (A,B,C,...) and female genotype
  counts must be supplied as a triangular matrix. See the examples below.
}
\value{
   \item{pval}{p-value}
   \item{tab}{table with allele counts}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{fisher.test}}}
\examples{
#
# bi-allelic autosomal
#

m <- c(AA=60,AB=96,BB=44) 
f <- c(AA=44,AB=97,BB=59) 
EAFtest <- EAFExact(m,f)

#
# bi-allelic X-chromosomal
#

males   <- c(A=392, B=212)
females <- c(AA=275, AB=296, BB=80)
EAFtest <- EAFExact(males,females,verbose=TRUE)

#
# tri-allelic autosomal
#

males   <- c(AA=20,AB=52,AC=34,BB=17,BC=51,CC=26)
females <- c(AA=28,AB=55,AC=33,BB=18,BC=50,CC=16)
EAFtest <- EAFExact(males,females,verbose=TRUE)

#
# tri-allelic X-chromosomal
#

males   <- c(A=15,B=17,C=24) 
females <- toTriangular(c(AA=4,AB=2,AC=13,BB=6,BC=19,CC=4)) 
EAFtest <- EAFExact(males,females,verbose=TRUE)

}
\keyword{htest}
