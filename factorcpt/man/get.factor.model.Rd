\name{get.factor.model}
\alias{get.factor.model}

\title{Factor model estimation via Principal Component Analysis}

\description{
Estimates the components of the factor structure
for an input time series, such as loadings and factors,
as well as estimating the number of factors.
}

\usage{
get.factor.model(x, max.q = NULL, q = NULL, bn = TRUE,
              bn.op = 2, normalisation = TRUE)
}

\arguments{
  \item{x}{
  input time series matrix, with each row representing a time series
  }
  \item{max.q}{
  see \code{max.q} in \code{\link{factor.seg.alg}}
}
  \item{q}{
  the number of factors; if \code{bn=FALSE}, \code{q} needs to be given
  }
  \item{bn}{
  if \code{bn=TRUE}, the factor number is estimated using the information
  criterion-based estimator of Bai and Ng (2002)
  }
  \item{bn.op}{
  an index number for the information criterion-based estimator of Bai and Ng (2002);
  the default value \code{bn.op=2} is used in Barigozzi, Cho & Fryzlewicz (2016)
}
  \item{normalisation}{
  if \code{normalisation=TRUE}, each row \code{x} is standardised prior to PCA
  }
}

\value{
  \item{lam}{an \code{n}-by-\code{(min(n, T)-1)} matrix containing the estimated loadings
  }
  \item{f}{a \code{(min(n, T)-1)}-by-\code{T} matrix containing the estimated factors
  }
  \item{norm.x}{
  if \code{normalisation=TRUE}, row-wise standardised \code{x}; else \code{norm.x=x}
  }
  \item{q.hat}{
  estimated number of factors
  }
  \item{max.q}{
  the maximum factor number used for factor number estimation
  }
  \item{ic}{
  information criterion values computed at a range of factor numbers from \code{0} to \code{max.q}
  }
}

\references{
J. Bai and S. Ng (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. \cr

M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint.
}

\author{
Haeran Cho
}

\examples{
n <- 50; T <- 200
e <- matrix(rnorm(n*T), nrow=n) # idiosyncratic components
r <- 3 # factor number
Lam <- matrix(rnorm(n*r, 1, 1), nrow=n) # loadings
f <- matrix(rnorm(r*T), nrow=r) # factors
chi <- e*0 # common component
chp <- T/2 # change-point
chi <- Lam\%*\%f
x <- chi + sqrt(r)*e

gfm <- get.factor.model(x)
gfm$q.hat
}
