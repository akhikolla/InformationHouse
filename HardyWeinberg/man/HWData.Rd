\name{HWData}
\alias{HWData}
\title{Generate genetic marker data in or out of Hardy-Weinberg Equilibrium}
\description{
   HWData generates samples of genotypic counts under various
   schemes. It mainly uses sampling from the multinomial distribution given
   Hardy-Weinberg allele frequencies.
}
\usage{
HWData(nm = 100, n = rep(100, nm), f = rep(0, nm), p = NULL, conditional
= FALSE, exactequilibrium = FALSE, pdist = "runif", x.linked = FALSE, nA
= NULL,
n.males=rep(round(0.5*n),nm), ...) 
}
\arguments{
  \item{nm}{The number of bi-allelic markers.}
  \item{n}{The sample sizes.}    
  \item{f}{The inbreeding coefficients (only for autosomal markers)}
  \item{p}{a vector of allele frequencies}
  \item{conditional}{if \code{TRUE} Haldane's distribution is used for sampling,
    if \code{FALSE} a multinomial distribution is used. Replaces
    parameter \code{pfixed} from the previous version of the package}
  \item{exactequilibrium}{generates data in exact HWE if set to
    \code{TRUE}}
  \item{pdist}{Take a random allele frequency from a uniform or beta
    distribution of \code{pfixed = FALSE} and \code{p} is not given.}
  \item{x.linked}{Simulated autosomal markers (\code{x.linked=FALSE},
    the default) or X-chromosomal markers (\code{x.linked=TRUE})}
  \item{nA}{A vector of minor allele counts, one for each marker. If not
    specified, it will be calculated from \code{p}}
  \item{n.males}{The number of males (only relevant if \code{x.linked = TRUE})}
  \item{\dots}{specific parameters for the uniform or beta}  
}
\details{
  The \code{exactequilibrium} option only takes effects for autosomal
  markers (\code{x.linked=FALSE}) and multinomial sampling
  (\code{conditional=FALSE}).

  Option \code{pfixed} is deprecated and replaced by \code{conditional}

  \code{HWData} returns a matrix of genotype counts, \code{nm} by 3 for
  autsomal markers or \code{nm} by 5 for X-chromosomal markers. Output
  is no longer supplied in the compositional form. Function
  \code{HWClo} can be used to convert the genotype counts to a composition.

  If the inbreeding coefficient is specified (\code{f}) it will only
  take effect for autosomal markers (\code{x.linked=FALSE}) and
  multinomial sampling (\code{conditional=FALSE}).
}
\value{
  \item{X}{A matrix containing the genotype counts.}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{HWTernaryPlot}}}
\examples{
nm <- 100
n <- 100
out <- HWData(nm,n)
}
\keyword{datagen}
