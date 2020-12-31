\name{roughness_K}
\alias{roughness_K}
\title{Kernel roughness}
\description{
  Computes the squared L2 norm, also known as the roughness, of a kernel implemented in FKSUM based on its beta coefficients. NB: the
  coefficients will first be normalised so that the kernel represents a density function. Equivalent to norm_K()^2
}
\usage{
roughness_K(beta)
}
\arguments{
  \item{beta}{numeric vector of positive coefficients.}
}
\value{
  positive numeric representing the squared L2 norm, or roughness of the kernel with coefficients beta/norm_const(beta).
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\keyword{file}
