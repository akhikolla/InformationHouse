\name{HWLRtest}
\alias{HWLRtest}
\title{Perform likelihood ratio test comparing two nested scenarios for a
  bi-allelic genetic variant, distinguishing the two sexes.
}
\description{
  Program \code{HWLRtest} performs a likelihood ratio test comparing two
  scenarios for an autosomal bi-allelic genetic variant. The scenarios
  concern Hardy-Weinberg proportions (HWP) and equality of allele
  frequencies (EAF) in both sexes. The different scenarios are described
  by Graffelman & Weir (2017). 
}
\usage{
HWLRtest(x, y, scene.null = "S1", scene.alt = "S6", verbose = TRUE, tracing = 0)
}
\arguments{
  \item{x}{Male genotype counts}
  \item{y}{Female genotype counts}
  \item{scene.null}{Scenario under the null hypothesis (E.g. "S1")}
  \item{scene.alt}{Scenario under the alternative hypothesis (E.g. "S6")}
  \item{verbose}{print output or not}
  \item{tracing}{Show tracing of the numeric likelihood maximization (1)
    or not (0).}
}
\details{
  The different scenarios are indicated with S1, S2, S3, S4, S6 and
  S6. S1 refers to Hardy-Weinber proportions and equality of allele
  frequencies. S2 refers to equality of allele frequencies and equality
  of inbreeding coefficients for the two sexes. S3 refers to equality of
  allele frequencies irrespective of HWP. S4 refers to HWP irrespective
  of allele frequencies. S5 refers to equality of inbreeding
  coefficients irrespective of allele frequencies. S6 is unrestricted.  
}
\value{
  \item{G2}{Likelihood ratio statistic}
  \item{df}{Degrees of freedom of the likelihood ratio statistic}
  \item{pval}{p-value}
}
\references{
 Graffelman, J. and Weir, B.S. (2017) On the testing of Hardy-Weinberg
 proportions and equality of allele frequencies in males and females at
 bi-allelic genetic
 markers. \href{http://dx.doi.org/10.1002/gepi.22079}{Genetic Epidemiology} pp. 1--15.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu}} 
\seealso{
\code{\link{HWAIC}}
}
\examples{
males <- c(AA=11,AB=32,BB=13) 
females <- c(AA=14,AB=23,BB=11)

#
# test EAF
#

lr1.out <- HWLRtest(males,females,scene.null="S3",scene.alt="S6")

#
# test EIC given EAF
#


lr2.out <- HWLRtest(males,females,scene.null="S2",scene.alt="S3")

#
# test HWP versus EIC, given EAF.
#

lr3.out <- HWLRtest(males,females,scene.null="S1",scene.alt="S2")

}
\keyword{htest}

