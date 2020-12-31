\name{HWABO}
\alias{HWABO}
\title{
Estimate allele frequencies and test for Hardy-Weinberg equilibrium with
a tri-allelic ABO system.
}
\description{
Function \code{HWABO} takes four genotype counts ("A","B","AB","OO") and
estimates the three allele frequencies using the EM algorithm.
}
\usage{
HWABO(x, p = c(1/3, 1/3, 1/3), maxiter = 50, tol = 1e-10, verbose = TRUE) 
}
\arguments{
  \item{x}{ a vector with genotype counts ("A","B","AB","OO").}
  \item{p}{ a vector with initial allele frequencies (by default (1/3,1/3,1/3)).}
  \item{maxiter}{ maximum number of iterations.}
  \item{tol}{ tolerance for convergence, 1e-10 by default}
  \item{verbose}{ print iteration history or not.}
}
\value{
    \item{pn }{ vector with estimated allele frequencies.}
    \item{It.hist }{ iteration history with log-likelihood.}
    \item{expected}{ expected genotype counts under HWE.}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
\code{\link{af}}
}
\examples{
   x <- c(fA=182,fB=60,nAB=17,nOO=176)
   al.fre <- HWABO(x)
#al2 <- HWABO(x,p=c(0.99,0.01,0.01),maxiter=25)
#al3 <- HWABO(x,p=c(0.01,0.99,0.01),maxiter=25)
#al4 <- HWABO(x,p=c(0.01,0.01,0.99),maxiter=25)
}
\keyword{misc}

