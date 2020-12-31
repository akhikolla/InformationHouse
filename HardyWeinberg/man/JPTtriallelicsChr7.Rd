\name{JPTtriallelicsChr7}
\alias{JPTtriallelicsChr7}
\docType{data}
\title{
Tri-allelic variants on chromosome 7 of the Japanese (JPT) sample of the 1000 genomes project 
}
\description{
\code{JPTtriallelics} contains six selected tri-allelic variants on 
chromosome 7 from the Japanese sample of the 1000 genomes project.
}
\usage{data("JPTtriallelicsChr7")}
\format{
  A data frame with 6 observations on the following 14 variables.
  \describe{
    \item{\code{id}}{RS identifier}
    \item{\code{pos}}{position in base pairs}
    \item{\code{mAA}}{number of AA males}
    \item{\code{mAB}}{number of AB males}
    \item{\code{mAC}}{number of AC males}
    \item{\code{mBB}}{number of BB males}
    \item{\code{mBC}}{number of BC males}
    \item{\code{mCC}}{number of CC males}
    \item{\code{fAA}}{number of AA females}
    \item{\code{fAB}}{number of AB females}
    \item{\code{fAC}}{number of AC females}
    \item{\code{fBB}}{number of BB females}
    \item{\code{fBC}}{number of BC females}
    \item{\code{fCC}}{number of CC females}
  }
}
\source{
  The \href{https://www.internationalgenome.org}{The 1000 genomes project}.
}
\references{
Graffelman, J. and Weir, B.S. (2017) Multi-allelic exact tests for Hardy-Weinberg equilibrium that
account for gender. doi: 10.1101/172874. Table 9.
}
\examples{
data(JPTtriallelicsChr7)
str(JPTtriallelicsChr7)
}
\keyword{datasets}
