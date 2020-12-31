\name{HWLratio}
\alias{HWLratio}
\title{Likelihood ratio test for Hardy Weinberg equilibrium}
\description{
  \code{HWLratio} performs the Likelihood ratio test for Hardy Weinberg
  equilibrium, both for autosomal and X-chromosomal markers.
}
\usage{
HWLratio(X, verbose = TRUE, x.linked = FALSE) 
}
\arguments{
  \item{X}{\code{X} a vector containing the genotypic counts (AA,AB,BB).}
  \item{verbose}{\code{verbose = TRUE} prints results, \code{verbose = FALSE}
    is silent.}
  \item{x.linked}{ \code{x.linked=FALSE} indicates the marker is
    autosomal (default), and \code{x.linked=TRUE} indicates it resides
    on the X-chromosome.}
}
\value{
\code{HWLratio} returns a list with the components:
  \item{Lambda }{the likelihood ratio}
  \item{G2 }{-2*log(Lambda)}
  \item{pval}{the p-value}
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer Associates, Massachusetts. See Chapter 3.}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWLratio(x,verbose=TRUE)
#
# Test for an X-chromsomal SNP.
#
rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
HW.test <- HWLratio(rs5968922,x.linked=TRUE,verbose=TRUE)
}
\keyword{htest}
