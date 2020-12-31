\name{post.cpts.analysis}
\alias{post.cpts.analysis}

\title{Display the results from change-point analysis and factor modelling}

\description{
Computes the minimum number of eigenvalues required so that
the given proportion of the variance of \code{x} over each segment,
defined by the change-points detected from the common components,
by varying the proportion and plots the number of eigenvalues.
}

\usage{
post.cpts.analysis(x, est.cps, cutoff.seq = seq(0.5, 0.95, by = 0.05), do.plot = TRUE)
}

\arguments{
  \item{x}{
  input time series
  }
  \item{est.cps}{
  a vector containing the change-points estimated for the common components
  }
  \item{cutoff.seq}{
  a sequence containing the proportions between zero and one
  }
  \item{do.plot}{
  if \code{do.plot=TRUE}, produces a plot of minimum numbers of eigenvalues required
  for accounting for the given proportion of variance, over time (\code{x}-axis)
  and the level of proportion (\code{y}-axis).
  }
}

\value{
a matrix containing the minimum numbers of eigenvalues for varying levels of proportions
over the segments defined by two neighbouring change-points in the common components
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint.
}
\author{Haeran Cho}

\examples{
n <- 50; T <- 200
e <- matrix(rnorm(n*T), nrow=n) # idiosyncratic components
r <- 3 # factor number
Lam <- matrix(rnorm(n*r, 1, 1), nrow=n) # loadings
f <- matrix(rnorm(r*T), nrow=r) # factors
chi <- e*0 # common component
chp <- T/2 # change-point
chi[, 1:chp] <- Lam\%*\%f[, 1:chp]
Lam <- Lam + matrix(rnorm(n*r, 0, sqrt(2)), nrow=n) # new loadings
chi[, (chp+1):T] <- Lam\%*\%f[, (chp+1):T]
x <- chi + sqrt(r)*e

post.cpts.analysis(x, 100)

}
