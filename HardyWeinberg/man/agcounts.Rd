\name{agcounts}
\alias{agcounts}
\title{
  Calculate allele and genotype counts for X-chromosomal markers
}
\description{
  Function \code{agcounts} determines sample size, minor are major allele counts,
  allele counts in females, numbers of males and females and allele
  frequencies for a vector of genotypes counts of an X-chromosomal markers.
}
\usage{
agcounts(x, verbose = FALSE)
}
\arguments{
  \item{x}{a vector of X-chromosomal genotype counts (A,B,AA,AB,BB)}
  \item{verbose}{print the counts if (\code{verbose = TRUE})}
}
\value{
  \item{n}{sample size}
  \item{nA}{number of A alleles}
  \item{nB}{number of B alleles}
  \item{nf}{number of females}
  \item{nm}{number of males}
  \item{nAf}{number of A alleles in females}
  \item{nBf}{number of B alleles in females}
  \item{nt}{total number of alleles}
  \item{fAA}{number of AA females}
  \item{fAB}{number of AB females}
  \item{fBB}{number of BB females}
  \item{pA}{overall A allele frequency}
  \item{pB}{overall B allele frequency}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{mac}},\code{link{maf}}}
\examples{
rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
counts <- agcounts(rs5968922)
}
\keyword{misc}

