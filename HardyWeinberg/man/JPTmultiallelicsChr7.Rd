\name{JPTmultiallelicsChr7}
\alias{JPTmultiallelicsChr7}
\docType{data}
\title{
  Multi-allelic autosomal variants of the Japanese population of the
  1000 genomes project
}
\description{
\code{JPTtriallelicsChrX} contains three selected multi-allelic variants
on chromosome 7 from the Japanese sample of the 1000 genomes project.
}
\usage{data("JPTmultiallelicsChr7")}
\format{List object with fields m4,f4; m5,f5; m6,f6;}
\details{
  The list object contains male and female genotype counts for 3
  multi-allelic variants on chromosome 7 of the JPT sample of the 1000 genomes project.
}
\source{
  The \href{https://www.internationalgenome.org}{The 1000 genomes project}.
}
\references{
Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that
account for gender. doi: 10.1101/172874. Table 10.
}
\examples{
data(JPTmultiallelicsChr7)
str(JPTmultiallelicsChr7)
}
\keyword{datasets}
