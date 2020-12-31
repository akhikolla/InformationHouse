\name{AFEItest}
\alias{AFEItest}
\title{ k-Sample AFEI Test of Equal Distributions}
\description{
Performs the distribution free exact k-sample test for equality of multivariate distributions in the HDLSS regime. This an aggregate test of the two sample versions of the FEI test over \eqn{\frac{k(k-1)}{2}} numbers of two-sample comparisons, and the test statistic is the minimum of these two sample FEI test statistics. Holm's step-down-procedure (1979) is applied for multiple testing.
}
\usage{
AFEItest(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", 
s_fn = 1, lb = 1, n_sts = 1000, alpha = 0.05)
}
\arguments{
\item{M}{\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes}
  \item{sizes}{vector of sample sizes}
  \item{randomization}{logical; if TRUE (default), randomization test and FALSE, non-randomization test}
  \item{clust_alg}{\code{"knwClustNo"}(default) or \code{"estclustNo"}; modified K-means algorithm used for clustering}
  \item{s_fn}{function required for clustering, 1 (default) for \eqn{1-\exp(-t)}, 2 for \eqn{\log(1+t)}, 3 for \eqn{t}}
  \item{lb}{each observation is partitioned into some numbers of smaller vectors of same length \eqn{lb}, default: \eqn{1}  }
  \item{n_sts}{number of simulation of the test statistic, default: \eqn{1000}}
  \item{alpha}{numeric, confidence level \eqn{\alpha}, default: \eqn{0.05}}
}
\value{
%%  ~Describe the value returned
   AFEItest returns a list containing the following items:
  \item{AFEIStat}{value of the observed test statistic}
  \item{AFCutoff}{cut-off of the test}
  \item{randomGamma}{randomized coefficient of the test}
  \item{decisionAFEI}{if returns \eqn{1}, reject the null hypothesis and if returns \eqn{0}, fails to reject the null hypothesis}
  \item{multipleTest}{indicates where two populations are different}
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020).  Distribution-free Exact k-sample Tests for High Dimensional Low Sample Size data Based on Clustering.

Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.

Sture Holm (1979). A simple sequentially rejective multiple test procedure, \emph{Scandinavian journal of statistics}, 65-70, doi:10.2307/4615733.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
  # muiltivariate normal distribution:
  # generate data with dimension d = 500
  set.seed(151)
  n1=n2=n3=n4=10
  d = 500
  I1 <- matrix(rnorm(n1*d,mean=0,sd=1),n1,d)
  I2 <- matrix(rnorm(n2*d,mean=0.5,sd=1),n2,d) 
  I3 <- matrix(rnorm(n3*d,mean=1,sd=1),n3,d) 
  I4 <- matrix(rnorm(n4*d,mean=1.5,sd=1),n4,d) 
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  #AFEI test:
  results <- AFEItest(M=X, sizes = c(n1,n2,n3,n4))
  
   ## outputs:
   results$AFEIStat
   #[1] 5.412544e-06

   results$AFCutoff
   #[1] 0.0109604

   results$randomGamma
   #[1] 0

   results$decisionAFEI
   #[1] 1
   
   results$multipleTest
   #  Population.1 Population.2 pvalues criticalValues
   #1            1            2       0    0.008333333
   #2            1            3       0    0.010000000
   #3            1            4       0    0.012500000
   #4            2            3       0    0.016666667
   #5            2            4       0    0.025000000
   #6            3            4       0    0.050000000

}
\keyword{ AFEItest function}