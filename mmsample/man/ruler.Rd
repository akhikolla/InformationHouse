\name{ruler}
\alias{ruler}
\title{Rcpp distance calculator}
\usage{
ruler(vR, uR, ciR)
}
\arguments{
\item{uR}{a vector of length k containing a list of values for all features (k) for the target. Numeric and dense.}
\item{vR}{an n x k matrix containing a matrix of values for all features (k) for all candidates (n). Numeric and dense.}
\item{ciR}{a square k x k matrix containing the inverted covariance matrix for all features (k). Numeric and dense.}
}
\description{
Returns a vector of distances from all rows in vR to the single row uR using ciR as the inverted covariance matrix.
}
\examples{
set.seed(123)
df <- data.frame(x = rpois(10, 20), y = rnorm(10, 50, 10))
cov_inv <- MASS::ginv(cov(df))
mmsample::ruler(as.matrix(df[2:10, ]), as.numeric(df[1, ]), cov_inv)
}
