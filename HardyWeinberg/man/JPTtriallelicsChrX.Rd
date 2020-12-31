\name{JPTtriallelicsChrX}
\alias{JPTtriallelicsChrX}
\docType{data}
\title{
Tri-allelic variants on the X-chromosome of the Japanese (JPT) sample of the 1000 genomes project 
}
\description{
\code{JPTtriallelicsChrX} contains five selected tri-allelic variants on the
X chromosome from the Japanese sample of the 1000 genomes project.
}
\usage{data("JPTtriallelicsChrX")}
\format{
  A data frame with 5 observations on the following 12 variables.
  \describe{
    \item{\code{id}}{Identifier of the polymorphism}
    \item{\code{pos}}{Position of the polymorphism in base pairs}
    \item{\code{chr}}{Chromosome}
    \item{\code{A}}{Number of males with A genotype}
    \item{\code{B}}{Number of males with B genotype}
    \item{\code{C}}{Number of males with C genotype}
    \item{\code{AA}}{Number of AA females}
    \item{\code{AB}}{Number of AB females}
    \item{\code{AC}}{Number of AC females}
    \item{\code{BB}}{Number of BB females}
    \item{\code{BC}}{Number of BC females}
    \item{\code{CC}}{Number of CC females}
  }
}
\source{
  The \href{https://www.internationalgenome.org}{The 1000 genomes project}.
}
\references{
Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that
account for gender. doi: 10.1101/172874. Table 6.
}
\examples{
data(JPTtriallelicsChrX)
str(JPTtriallelicsChrX)
}
\keyword{datasets}
