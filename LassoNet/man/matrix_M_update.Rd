\name{matrix.M.update}
\alias{matrix.M.update}
\title{
Updates connection sign matrix.
}
\description{
Updates M using relation \eqn{(M)_{ij} = -\hat{\xi}_{ij} |(M_1)|_{ij}}.
}
\usage{matrix.M.update(M, xi)}
\arguments{
\item{M}{penalty matrix}
\item{xi}{estimated \eqn{\hat{\xi}_{ij}} matrix}
}
\details{Updates M}
\value{
\item{M}{updated M}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{Weber, M., Striaukas, J., Schumacher, M., Binder, H. "Network-Constrained Covariate Coefficient and Connection Sign Estimation" (2018) <doi:10.2139/ssrn.3211163>}

\examples{
p<-100
M<-diag(p)
xi<-matrix(rnorm(p*p), p, p)
matrix.M.update(M,xi)
}
