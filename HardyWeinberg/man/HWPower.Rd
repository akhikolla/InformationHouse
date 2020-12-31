\name{HWPower}
\alias{HWPower}
\title{
  Compute the power of a test for Hardy-Weinberg equilibrium.
}
\description{
\code{HWPower} is a function that computes the power of a test for
Hardy-Weinberg equilibrium.
}
\usage{
HWPower(n = 100, nA = 100, pA = 0.5, y = c(AA=25,AB=50,BB=25),
alpha = 0.05, theta = 4, f = NULL, test = "exact",
alternative = "two.sided", pvaluetype = "selome", cc = 0.5)
}
\arguments{
  \item{n}{The sample size}
  \item{nA}{The minor allele count}
  \item{pA}{The minor allele frequency}
  \item{y}{A sample of genotype counts (AA,AB,BB)}
  \item{alpha}{The significance level (0.05 by default)}
  \item{theta}{The degree of disequilibrium (\code{theta} = 4 is
    equilibrium, \code{theta} > 4 is heterozygote excess, \code{theta} <
    4 is heterozygote dearth)}
  \item{f}{The inbreeding coefficient. Overrules \code{theta} if specified.}
  \item{test}{The type of test for which power is to be computed. Can be
    "exact" (default) or "chisq" (chi-square)}
  \item{alternative}{The nature of the alternative hypothesis
    ("two.sided" (default), "greater" or "less")}
  \item{pvaluetype}{The type of p-value used in an exact test
    ("selome", "dost" or "midp")}
  \item{cc}{Continuity correction parameter for the chi-square test (0.5
    by default)}
}
\details{
\code{HWPower} uses the Levene-Haldane distribution (distribution of the
number of heterzygotes given the minor allele count) for computing
power.

\code{HWPower} can be used in three different way. In principle, the
power is calcuted on the basis of the sample size (\code{n}) and the minor
allele count (\code{nA}). Alternatively, the user may specifiy sample
size (\code{n}) and minor allele frequency (\code{pA}). Finally, power
can also be calculated directly from a sample of genotype counts. In
that case the calculated power is the power for a sample of the given
sample size and minor allele count. The three ways to use \code{HWPower}
are illustrated in the example section. 

}
\value{
  if \code{test} = "exact" the power of the exact test is computed for
  the given significance level and minor allele count.

  if \code{test} = "chisq" the power of the chi-square test is computed for
  the given significance level and minor allele count.
}
\references{
Graffelman, J. and Moreno, V. (2013) The Mid p-value in exact tests for
Hardy-Weinberg proportions. \emph{Statistical Applications in Genetics
  and Molecular Biology} 12(4): 433-448. 

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3):
1-23. \url{https://www.jstatsoft.org/v64/i03/}.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{HWExact}}}
\examples{
pw.chisq <- HWPower(n=100,nA=100,alpha=0.05,test="chisq",theta=16)
print(pw.chisq)
pw.exact <- HWPower(n=100,nA=100,alpha=0.05,test="exact",theta=16,pvaluetype="selome")
print(pw.exact)
pw.exact <- HWPower(n=100,nA=100)
print(pw.exact)
pw.exact <- HWPower(n=100,pA=0.5)
print(pw.exact)
pw.exact <- HWPower(y=c(AA=25,AB=50,BB=25))
print(pw.exact)
}
\keyword{htest}
