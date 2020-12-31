\name{HWExact}
\alias{HWExact}
\title{ Exact test for Hardy-Weinberg equilibrium }
\description{
  \code{HWExact} performs an exact test for Hardy-Weinberg equilibrium
}
\usage{
HWExact(X, alternative = "two.sided", pvaluetype = "selome", eps=1e-10, x.linked =
FALSE, verbose = TRUE)
}
\arguments{
  \item{X}{ vector with the genotype counts AA, AB, BB }
  \item{alternative}{\code{two.sided} (default) will perform a two-sided
    test where both an excess and a dearth of heterozygotes count as
    evidence against HWE. \code{less} is a one-sided test where only
    dearth of heterozygotes counts a evidence against HWE,
    \code{greater} is a one-sided test where only excess of
    heterozygotes counts as evidence against HWE.}
  \item{pvaluetype}{if \code{pvaluetype} is set to \code{dost}
then the p-value of a two-sided test is computed as twice the tail area
of a one-sided test. When set to \code{selome}, the p-value
is computed as the sum of the probabilities of all samples less or
equally likely as the current sample. When set to \code{midp}, the
p-value is computed as half the probability of the current sample +
the probabilities of all samples that are more extreme.}
  \item{x.linked}{ \code{x.linked=FALSE} indicates the marker is
    autosomal (default), and \code{x.linked=TRUE} indicates it resides
    on the X-chromosome.}
  \item{eps}{a tolerance that can be set for comparing probabilities in order to include
    tied outcomes}
   \item{verbose}{ print results or not.}
}
\details{
   \code{HWExact} uses the recursion equations described by Wigginton
   et. al.

   For testing large sets of bi-allelic variants, use the faster code in \code{\link{HWExactStats}}.
   
   For large samples, \code{HWExact} may give the error message:
   "evaluation nested too deeply: infinite recursion". This can usually
   be resolved by increasing R's limit on nested expressions with
   \code{options(expressions=10000)} or a higher limit. With higher
   limits, the error message "protect(): protection stack overflow" can
   occur. This error can usually be resolved by increasing R's
   protection stack with the command line option
   \code{--max-ppsize 100000} or higer values. However, with such large
   samples the exact test will give virtually the same result as a
   chi-square test, and it may be easier to use \code{HWChisq} in these
   circumstances.
}
\value{
  \item{pval}{ p-value of the exact test }
  \item{prob}{ probabilities of all possible samples with the same
   sample size and minor allele count}
  \item{pofthesample}{ probability of the observed sample }
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.
  
  Wigginton, J.E., Cutler, D.J. and Abecasis, G.R. (2005) A note on exact
  tests of Hardy-Weinberg equilibrium, American Journal of Human
  Genetics (76) pp. 887-893.

  Graffelman, J. and Moreno, V. (2013) The mid p-value in exact tests
  for Hardy-Weinberg equilibrium, Statistical Applications in Genetics
  and Molecular Biology 12(4) pp. 433-448.}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisq}}, \code{\link{HWExactStats}} }
\examples{
#
# Example for an autosomal marker using the standard exact p-value
#
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWExact(x,verbose=TRUE)
#
# Example for an autosomal marker using the mid p-value
#
HW.test <- HWExact(x,verbose=TRUE,pvaluetype="midp")
#
# Example x-linked marker
#
rs5968922  <- c(A=392, B=212, AA=275, AB=296, BB=80 )
HWExact(rs5968922,x.linked=TRUE,verbose=TRUE) 
}
\keyword{htest}
