\name{gMADD_DI}
\alias{gMADD_DI}
\title{Modified K-Means Algorithm by Using a New Dissimilarity Measure, MADD and DUNN Index
}
\description{Performs modified K-means algorithm by using a new dissimilarity measure, called MADD and DUNN index, and provides estimated cluster (class) labels or memberships and corresponding DUNN index of the observations. 
}
\usage{
gMADD_DI(s_fn, kmax, lb, M)
}
\arguments{
  \item{s_fn}{function required for clustering, 1 for \eqn{1-\exp(-t)}, 2 for \eqn{\log(1+t)}, 3 for \eqn{t}
}
  \item{kmax}{ maximum value of total number of clusters to estimate total number of the classes in the whole observations 
}
  \item{lb}{each observation is partitioned into some numbers of smaller vectors of same length \eqn{lb}
}
  \item{M}{\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes
}
}
\details{
DUNN index is used for cluster validation, but here we use it to estimate total number of cluster \eqn{k} by \eqn{\hat k = argmax_{2\le k' \le k^*}DI(k')}. Here \eqn{DI(k')} represents the DUNN index and we use \eqn{k^*=2*k}.
}
\value{
a \eqn{kmax \times (n+1)} matrix of the estimated cluster (class) labels and corresponding DUNN indexes of observations
}
\note{
The result of this gMADD_DI function is a matrix. The 1st row of this matrix doesn't provide anything about estimated class labels or DUNN index of observations since the DUNN index is only defined for \eqn{k\ge 2}. The last column of this matrix represents the DUNN indexes. The estimated cluster labels of observations are calculated by finding out the corresponding row of maximum DUNN index.
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020).  Distribution-free Exact k-sample Tests for High Dimensional Low Sample Size data Based on Clustering.

Soham Sarkar and Anil K Ghosh (2019). On perfect clustering of high dimension, low sample size data, \emph{IEEE transactions on pattern analysis and machine intelligence}, doi:10.1109/TPAMI.2019.2912599.

Joseph C Dunn (1973). A fuzzy relative of the isodata process and its use in detecting compact well-separated clusters, doi:10.1080/01969727308546046.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
  # Modified K-means algorithm:
  # muiltivariate normal distribution
  # generate data with dimension d = 500
  set.seed(151)
  n1=n2=n3=n4=10
  d = 500
  I1 <- matrix(rnorm(n1*d,mean=0,sd=1),n1,d)
  I2 <- matrix(rnorm(n2*d,mean=0.5,sd=1),n2,d) 
  I3 <- matrix(rnorm(n3*d,mean=1,sd=1),n3,d) 
  I4 <- matrix(rnorm(n4*d,mean=1.5,sd=1),n4,d) 
  n_cl <- 4
  N <- n1+n2+n3+n4
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  dvec_di_mat <-  gMADD_DI(1,2*n_cl,1,X)
  est_no_cl <- which.max(dvec_di_mat[ ,(N+1)])
  dvec_di_mat[est_no_cl,1:N]

   ## outputs:
   #[1] 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3

}
\keyword{ gMADD_DI function}