\name{var_K}
\alias{var_K}
\title{Variance of a kernel}
\description{
  Computes the variance of the random variable whose density is given by a kernel implemented in FKSUM, with coefficients beta.
  NB: coefficients will first be normalised so that the kernel is a density function.
}
\usage{
var_K(beta)
}
\arguments{
  \item{beta}{positive numeric vector of kernel coefficients.}
}
\value{
  A positive numeric value representing the variance of the random variable with density given by the kernel.
}
\keyword{file}
