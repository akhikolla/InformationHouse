\name{norm_K}
\alias{norm_K}
\title{The L2 norm of a kernel}
\description{
  Computes the L2 norm of a kernel implemented in FKSUM, based on its coefficients. NB: the coefficients will first be normalised so that
  the kernel represents a density function. Equivalent to sqrt(roughness_K())
}
\usage{
norm_K(beta)
}
\arguments{
  \item{beta}{numeric vector of positive coefficients.}
}
\value{
  positive numeric representing the L2 norm of the kernel with coefficients beta/norm_const(beta).
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\keyword{file}
