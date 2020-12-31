\name{HWPosterior}
\alias{HWPosterior}
\title{
Calculation of posterior probabilities and Bayes factors for
Hardy-Weinberg tests at X-chromosomal variants.
}
\description{
Function \code{HWPosterior} calculates posterior probabilities and Bayes
factors for tests for Hardy-Weinberg equilibrium of autosomal and X-chromosomal
variants.
}
\usage{
HWPosterior(X, verbose = TRUE, prior.af = c(0.5,0.5), prior.gf =
            c(0.333,0.333,0.333), x.linked = FALSE, precision = 0.05) 
}
\arguments{
  \item{X}{A vector of genotype counts. The order c(A,B,AA,AB,BB) is
    assumed. Differently ordered vectors can be supplied but then elements must
    be labeled by their genotype}
  \item{verbose}{prints results if \code{verbose = TRUE}}
  \item{prior.af}{Beta prior parameters for male and female allele frequencies}
  \item{prior.gf}{Dirichlet prior parameters for female genotype
    frequencies}
  \item{x.linked}{logical indicating whether the variant is autosomal or
    X-chromosomal}
  \item{precision}{precision parameter for marginal likelihoods that
    require numeric integration}
}
\details{For X-chromosomal variants, four possible models are considered, and the posterior
  probabilities and Bayes factors for each model are calculated.

  For autosomal variants, ten possible scenarios are considered, and the
  posterior probabilities for all models are calculated.

  In general, default Dirichlet priors are used for genotype
  frequencies, and beta prior are used for allele frequencies.
}
\value{
  For X-chromosomal variants, a matrix with posterior probabilities and
  Bayes factors will be produced.
  For autosomal variants, a vector of posterior probabilities is produced.
}
\references{
  Puig, X., Ginebra, J. and Graffelman, J. (2017) A Bayesian test for
  Hardy-Weinberg equilibrium of bi-allelic X-chromosomal markers. To
  appear in Heredity.
}
\author{
  Xavi Puig \email{xavier.puig@upc.edu} and
  Jan Graffelman \email{jan.graffelman@upc.edu}
}
\seealso{ \code{\link{HWChisq}}, \code{\link{HWExact}}, \code{\link{HWExactStats}} }
\examples{
#
# An X-chromosomal example
#
x <- c(A=43,B=13,AA=26,AB=19,BB=3)
out <- HWPosterior(x,verbose=TRUE,x.linked=TRUE)
#
# An autosomal example
#
data(JPTsnps)
post.prob <- HWPosterior(JPTsnps[1,],x.linked=FALSE)
}
\keyword{htest}

