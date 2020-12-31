\name{recode}
\alias{recode}
\title{Recode genotype information}
\description{
   function \code{recode} recodes bi-allelic genetic marker information
   expressed as strings (e.g. "AA", "AB", "BB") into numerical form.
}
\usage{recode(X, alleles, values = c(0, 1, 2), pos1 = 1, pos2 = 3, minor
= FALSE, verbose = FALSE)}
\arguments{
  \item{X}{A matrix or dataframe of bi-allelic markers, individuals in
    rows, markers in columns}
  \item{alleles}{a vector with the alleles for each marker
    (e.g. c("A/T", "A/G", etc))}
  \item{values}{a vector of numerica values for AA, AB and BB, ((0,1,2)
  by default).}
  \item{pos1}{position of the first allele in the allele string (1 by default).}
  \item{pos2}{position of the second allele in the allele string (3 by
    default).}
  \item{minor}{coding is according to the number of copies of the minor
    allele. if \code{minor = TRUE}, the value of 2 reflects two copies
    of the minor allele, and the value 0 reflects no copies of the minor
    allele.}
  \item{verbose}{print progress on the conversion or not.}
}
\details{
  \code{recode} is written for bi-allelic marker data
  only. Heterozygotes may be coded both as AB or BA. By default, the
  second allele specified (e.g. "T" in "A/T") is counted in the
  recoding, and homozygotes AA are coded as 0 and homozygotes TT as 2. 
}
\value{A numerical matrix, individuals in rows, markers in columns}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{MakeCounts}}}
\examples{
SNP1 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP2 <- c("CG","GG","CC","GG","GG","CG","CG","CG","CG")
SNP3 <- c("AA","AA","AA","AG","AA","AG","AA","AA","AA")
SNP4 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP5 <- c("CC","CC","CC","CC","CC","CC","CT","CT","CT")
X <- cbind(SNP1,SNP2,SNP3,SNP4,SNP5)
Y <- recode(X,c("A/G","C/G","A/G","A/G","C/T"))
print(Y)
}
\keyword{misc}

