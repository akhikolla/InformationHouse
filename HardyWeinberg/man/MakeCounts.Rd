\name{MakeCounts}
\alias{MakeCounts}

\title{Create genotype counts from bi-allelic marker data}
\description{
\code{MakeCounts} creates a matrix of genotype counts, with one row for
each bi-allelic marker, containing 4 columns with the counts AA, AB, BB
and NA (missings) respectively
}
\usage{
MakeCounts(X, alleles, pos1 = 1, pos2 = 3, coding = c(AA=0,AB=1,BB=2), sep = "")
}
\arguments{
  \item{X}{A matrix or dataframe with bi-allelic genotyping information,
    markers in columns, individuals in rows}
  \item{alleles}{a vector of alleles for each marker
    (e.g. c("A/T","A/G",...)). Only relevanit if \code{X} is a matrix
    with text entries.}
  \item{pos1}{position of the first allele in the allele string (1 by default)}
  \item{pos2}{position of the second allele in the allele string (3 by
    default)}
  \item{coding}{indicates how homozygotes and heterozygote are coded as
    numbers. Only relevant if \code{X} is a matrix with numeric
    entries.}
  \item{sep}{allele separator character for genotype data in text format ("" for AA; "/" for "A/A")}
}
\details{
  \code{MakeCounts} is thought for bi-allelic marker data only. Missings
  are should be coded by NA. It produces the right input for
  \code{HWTernaryPlot}.

  Heterozygotes may be coded in the data as "AB" or "BA". Both entries
  will be counted as a heterozygote. 
}
\value{A matrix of 4 columns}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{HWTernaryPlot}}}
\examples{
SNP1 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP2 <- c("CG","GG","CC","GG","GG","CG","CG","CG","CG")
SNP3 <- c("AA","AA","AA","AG","AA","AG","AA","AA","AA")
SNP4 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP5 <- c("CC","CC","CC","CC","CC","CC","CT","CT","CT")
X <- cbind(SNP1,SNP2,SNP3,SNP4,SNP5)
Y <- MakeCounts(X,c("A/G","C/G","A/G","A/G","C/T"))
print(Y)
W <- matrix(sample(c(0,1,2,NA),100,replace=TRUE),ncol=5)
Z <- MakeCounts(W,coding=c(0,1,2))
}
\keyword{misc}

