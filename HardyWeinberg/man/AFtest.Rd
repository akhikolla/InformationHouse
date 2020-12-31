\name{AFtest}
\alias{AFtest}
\title{Exact test of equality of allele frequencies for males and females
}
\description{
Function \code{AFtest} tests equality of allele frequencies for males
and females for bi-allelic marker data by means of a Fisher exact test.  
}
\usage{
AFtest(x, verbose = TRUE, ...)
}
\arguments{
  \item{x}{a vector containg the genotypic counts c(A,B,AA,AB,BB) for a
    bi-allelic X-chromosomal markers.}
  \item{verbose}{verbose = TRUE prints results, verbose = FALSE is silent.}
  \item{\dots}{additional arguments for function \code{fisher.test}.}
}
\details{
  Function \code{AFtest} constructs the contingency table of sex by
  allele, and call \code{fisher.test} to test for equality of allele
  frequencies. The test assumes Hardy-Weinberg equilibrium.
}
\value{
  \item{AC}{Two-way table of sex by allele}
  \item{pval}{p-value of the test}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}}, \code{\link{HWExact}} }
\examples{
  rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
  AFtest(rs5968922)
}
\keyword{htest}

