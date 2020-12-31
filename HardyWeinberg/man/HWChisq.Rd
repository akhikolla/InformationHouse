\name{HWChisq}
\Rdversion{1.1}
\alias{HWChisq}
\title{Chi square tests for Hardy Weinberg equilibrium}
\description{
  \code{HWChisq} performs the chi-square test for Hardy-Weinberg
  equilibrium both for autosomal and X-chromosomal markers.  
}
\usage{
HWChisq(X, cc = 0.5, verbose = TRUE, x.linked = FALSE, phifixed = NULL) 
}
\arguments{
  \item{X}{\code{X} a vector containg the genotypic counts (AA,AB,BB for
    autosomal markers c(A,B,AA,AB,BB) for X-chromosomal markers).}
  \item{cc}{\code{cc} the continuity correction parameter (default
    \code{cc = 0.5}).}
  \item{verbose}{ \code{verbose = TRUE} prints results, \code{verbose = FALSE}
    is silent.}
  \item{x.linked}{ \code{x.linked=FALSE} indicates the marker is
    autosomal (default), and \code{x.linked=TRUE} indicates it resides
    on the X-chromosome.}
  \item{phifixed}{(For X-chromosomal markers only)
    \code{phifixed=NULL} indicates that the fraction of males
    (females) should be estimated from the data (default). If set to any
  other value (e.g. \code{phifixed=0.5}) then the sample is assumed to
  come from a population with the specified fraction of males.}
}
\details{
   \code{HWChisq} does a chi-square test for Hardy-Weinberg equilibrium,
   and by default applies a continuity correction. For extreme allele
   frequencies, the continuity correction can lead to excessive type 1
   error rates, and is better turned off in that case. The continuity
   correction can be turned off by specifying \code{cc=0}.

   \code{HWChisq} can do the chi-square test for both autosomal and
   X-chrosomal markers. By setting \code{x.linked = TRUE} the marker
   will be assumed to be on the X-chromosome, and the count vector
   supplied should have 5 elements instead of 3 elements for an
   autosomal marker. For X-chromsomal markers argument \code{phifixed}
   is in general best left to its default value (\code{NULL}). Only in
   specific situations where the theoretical population sex ratio is known (e.g. in
   simulation studies where a universe with known gender ratio is
   sampled) \code{phifixed} could be set to the theoretical ratio of interest.   

   When \code{alternative} is set to \code{less}, a one-sided test for
   against a negative inbreeding coefficient (heterozygote excess) is
   performed. When \code{alternative} is set to \code{greater} a one-sided test for
   against a positive inbreeding coefficient (lack of heterozygotes) is
   performed.
   
}
\value{
\code{HWChisq} returns a list with the components:
  \item{chisq }{value of the chi-square statistic. NA is returned if the marker is monomorphic.}
  \item{pval }{p-value of the chi-square test for Hardy-Weinberg equilibrium.}
  \item{D }{Half the deviation from Hardy-Weinberg equilibrium for the AB genotype.}
  \item{p }{the allele frequency of A.}
  \item{f }{the inbreeding coefficient.}
  \item{expected }{the expected counts under Hardy-Weinberg equilibrium} 
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.

  For the chi-square test for X-linked markers:

  Graffelman, J. & Weir, B.S. (2016) Testing for Hardy-Weinberg
  equilibrium at bi-allelic genetic markers on the X chromosome. 
  \emph{Heredity} 116(6) pp. 558--568. \url{http://dx.doi.org/10.1038/hdy.2016.20}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisqStats}} }
\examples{
#
# Test for an autosomal blood group marker
#
x <- c(MM=298,MN=489,NN=213)
HW.test <- HWChisq(x,verbose=TRUE)
#
# Same test without continuity correction
#
HW.test <- HWChisq(x,cc=0,verbose=TRUE)
#
# Test for an X-chromsomal SNP.
#
rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
HW.test <- HWChisq(rs5968922,cc=0,x.linked=TRUE,verbose=TRUE)
}
\keyword{htest}
