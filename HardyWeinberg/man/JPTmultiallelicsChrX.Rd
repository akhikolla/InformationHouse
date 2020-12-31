\name{JPTmultiallelicsChrX}
\alias{JPTmultiallelicsChrX}
\docType{data}
\title{
  Multi-allelic X-chromosomal variants of the Japanese population of the
  1000 genomes project
}
\description{
\code{JPTtriallelicsChrX} contains four selected multi-allelic variants on the
X chromosome from the Japanese sample of the 1000 genomes project.
}
\usage{data("JPTmultiallelicsChrX")}
\format{List object with fields m4,f4; m5,f5; m6,f6; m7,f7}
\details{
  The list object contains male and female genotype counts for four
  multi-allelic variants of the JPT sample of the 1000 genomes project.
}
\source{
  The \href{https://www.internationalgenome.org}{The 1000 genomes project}.
}
\references{
Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that
account for gender. doi: 10.1101/172874. Table 7.
}
\examples{
data(JPTmultiallelicsChrX)
m4 <- JPTmultiallelicsChrX$m4
f4 <- JPTmultiallelicsChrX$f4
}
\keyword{datasets}
