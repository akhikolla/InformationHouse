\name{TSIXTriAllelics}
\docType{data}
\alias{TSIXTriAllelics}
\title{
  Tri-allelic polymorphisms on the X chromosome of the TSI population
}
\description{
  This dataframe contains genotype counts for six three-allelic polymorphisms (A,B,C) on chromosome X of a sample of individuals from the TSI population (Tuscany, Italy) of the 1,000 genomes project.
}
\usage{data(TSIXTriAllelics)}
\format{
  A data frame with 6 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{Identifier of the polymorphism}
    \item{\code{Ma}}{Male A genotype count}
    \item{\code{Mb}}{Male B genotype count}
    \item{\code{Mc}}{Male C genotype count}
    \item{\code{Faa}}{Female AA genotype count}
    \item{\code{Fab}}{Female AB genotype count}
    \item{\code{Fac}}{Female AC genotype count}
    \item{\code{Fbb}}{Female BB genotype count}
    \item{\code{Fbc}}{Female BC genotype count}
    \item{\code{Fcc}}{Female CC genotype count}
  }
}
\source{
  Data taken from the 1,000 genomes project at www.internationalgenome.org
}
\references{
  Graffelman, J. and Ortoleva, L. (2020) A network algorithm for the X chromosomal exact test for Hardy-Weinberg equilibrium with multiple alleles. Under review.
}
\examples{
data(TSIXTriAllelics)
}
\keyword{datasets}
