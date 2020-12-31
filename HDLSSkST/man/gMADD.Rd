\name{gMADD}
\alias{gMADD}
\title{Modified K-Means Algorithm by Using a New Dissimilarity Measure, MADD
}
\description{Performs modified K-means algorithm by using a new dissimilarity measure, called MADD, and provides estimated cluster (class) labels or memberships of observations.
}
\usage{
gMADD(s_fn, n_clust, lb, M)
}
\arguments{
  \item{s_fn}{function required for clustering, 1 for \eqn{1-\exp(-t)}, 2 for \eqn{\log(1+t)}, 3 for \eqn{t}
}
  \item{n_clust}{total number of the classes in the whole observations
}
  \item{lb}{each observation is partitioned into some numbers of smaller vectors of same length \eqn{lb} 
}
  \item{M}{\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes
}
}

\value{
a vector of length n of estimated cluster (class) labels of observations
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020). Distribution-free Exact k-sample Tests for High Dimensional Low Sample Size data Based on Clustering.

Soham Sarkar and Anil K Ghosh (2019). On perfect clustering of high dimension, low sample size data, \emph{IEEE transactions on pattern analysis and machine intelligence}, doi:10.1109/TPAMI.2019.2912599.

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
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  gMADD(1,n_cl,1,X)
  
   ## outputs:
   #[1] 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3

}
\keyword{ gMADD function}