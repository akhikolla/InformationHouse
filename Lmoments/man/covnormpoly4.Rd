\name{covnormpoly4}
\alias{covnormpoly4}
\title{Covariance matrix of the parameters of the normal-polynomial quantile mixture}
\description{
Estimates covariance matrix of the four parameters of normal-polynomial quantile mixture}
\usage{
covnormpoly4(data)
}
\arguments{
  \item{data}{vector of observations}
}
\value{covariance matrix of the four parameters of normal-polynomial quantile mixture}

\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.
}
\seealso{\code{\link{Lmomcov}} for covariance matrix of L-moments, 
\code{\link{dnormpoly}} for the normal-polynomial quantile mixture and
\code{\link{data2normpoly4}} for the estimation of the normal-polynomial quantile mixture.}
\author{Juha Karvanen \email{juha.karvanen@iki.fi}}
\keyword{distribution}
\keyword{robust}


