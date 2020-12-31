\name{HWAIC}
\alias{HWAIC}
\title{Compute Akaike's Information Criterion (AIC) for HWP and EAF models}
\description{
  Function \code{HWAIC} calculates Akaike's Information Criterion for
  ten different models that describe a bi-allelic genetic variant: M11:
  Hardy-Weinberg proportions and equality of allele frequencies in the
  sexes (HWP & EAF); M12: EAF and HWP in males only; M13: EAF and HWP in
  females only; M14: EAF and equality of inbreeding coefficients in
  the sexes (EIC); M15: EAF only; M21: HWP in both sexes; M22: HWP for
  males only; M23: HWP for females only; M24: EIC only; M25: None of the previous.
}
\usage{
HWAIC(x, y, tracing = 0, tol = 0.000001)
}
\arguments{
  \item{x}{Male genotype counts (AA,AB,BB)}
  \item{y}{Female genotype counts (AA,AB,BB)}
  \item{tracing}{Activate tracing in the maximization of some
    likelihoods (0=no tracing; 1:tracing)}
  \item{tol}{tolerance for iterative maximization of some likelihoods}
}
\details{
  The log-likelihood for the six models is calculated. For two models (C
  and E) this is done numerically using package \code{RSolnp}.
}
\value{A named vector containing 6 values for AIC}
\references{
 Graffelman, J. and Weir, B.S. (2017) On the testing of Hardy-Weinberg
 proportions and equality of allele frequencies in males and females at
 bi-allelic genetic
 markers. \href{http://dx.doi.org/10.1002/gepi.22079}{Genetic Epidemiology} pp. 1--15.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu}} 
\seealso{\code{\link{HWLRtest}}}
\examples{
males <- c(AA=11,AB=32,BB=13) 
females <- c(AA=14,AB=23,BB=11) 
stats <- HWAIC(males,females)
print(stats)
}
\keyword{misc}

