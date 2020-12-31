\name{HWPerm}
\alias{HWPerm}
\title{
Permutation test for Hardy-Weinberg equilibrium
}
\description{
Function \code{HWPerm} does a permutation test for Hardy-Weinberg
equilibrium using a user-supplied test statistic.
}
\usage{
HWPerm(x, nperm = 17000, verbose = TRUE, x.linked = FALSE,
FUN = ifelse(x.linked,Chisquare.x,Chisquare), eps=1e-10, ...) 
}
\arguments{
  \item{x}{A vector of genotype counts (AA,AB,BB)}
  \item{nperm}{The number of permutations}
  \item{verbose}{\code{verbose} = TRUE will print results,
    \code{verbose} = FALSE is silent.}    
  \item{x.linked}{ \code{x.linked=FALSE} indicates the marker is
    autosomal (default), and \code{x.linked=TRUE} indicates it resides
    on the X-chromosome.}
  \item{FUN}{An function call for calculating the test statistic for HWE
    (see examples below)}
  \item{eps}{Tolerance for comparison of floating point numbers (1e-10 by default)}
  \item{\dots}{Additional parameters for the function call argument \code{FUN}}
}
\details{
The set of alleles for the observed sample is permuted. Consequently, the test is
conditional on allele frequency.
}
\value{
\code{HWPerm} returns a list with the components:
  \item{stat }{value of the chosen test statistic for the observed sample.}
  \item{pval }{p-value of the permutation test.}
}
\references{
Ziegler, A. & K\"{o}nig, I.R. (2006) A statistical approach to genetic
epidemiology. Wiley.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{HWChisq}},\code{\link{HWExact}},\code{\link{HWLratio}}}
\examples{
x <- c(MM=298,MN=489,NN=213)
\dontrun{
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE) # uses default chi-square statistic
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE,function(z)
HWChisq(z)$chisq,cc=0.5) # uses chi-square statistic with continuity correction.
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE,function(y) HWLratio(y)$G2)
# uses likelihood ratio statistic.
HWPerm(x,nperm=10000,verbose=TRUE,function(y) 1-HWExact(y)$pval)
# uses exact test p-value
#
# Permutation test for a marker on the X chromosome
#
rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
HW.test <- HWPerm(rs5968922,nperm=10000,x.linked=TRUE,verbose=TRUE)
}
}
\keyword{htest}

