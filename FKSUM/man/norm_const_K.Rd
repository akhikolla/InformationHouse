\name{norm_const_K}
\alias{norm_const_K}
\title{Normalising constant for kernels in FKSUM}
\description{
  Computes the normalising constant for kernels implemented in the FKSUM package, to convert kernel with arbitrary beta coefficients to one which is a probability density. That is the kernel with coefficients equal to beta/norm_const_beta(beta) has unit integral
}
\usage{
norm_const_K(beta)
}
\arguments{
  \item{beta}{numeric vector of positive coefficients.}
}
\value{
  positive numeric normalising constant
}
\references{
Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\keyword{file}
