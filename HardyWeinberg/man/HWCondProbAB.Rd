\name{HWCondProbAB}
\alias{HWCondProbAB}
\title{ Compute probability of a genotypic sample  }
\description{
  Computes the probability of a particular genotypic sample given the
  allele count, sample size and number of heterozygotes.
}
\usage{
HWCondProbAB(n, nA, nAB)
}
\arguments{
  \item{n}{ \code{n} is the total sample size (total number of individuals)}
  \item{nA}{ \code{nA} is the number of A alleles in the sample }
  \item{nAB}{ \code{nAB} is the number of heterozygotes in the sample }
}
\value{
  \item{p}{probability of the particular sample}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{HWExact}}}
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
n <- sum(x)
nM <- 2*x[1]+x[2]
nMN <- x[2]
p <- HWCondProbAB(n,nM,nMN)
}
\keyword{misc}

