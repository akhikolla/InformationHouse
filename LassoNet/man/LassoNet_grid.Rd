\name{lasso.net.grid}
\alias{lasso.net.grid}
\title{
Estimates coefficients and connection signs over the grid of values  of penalty parameters \eqn{\lambda}1 and \eqn{\lambda}2.
}
\description{
Fits network regressions over the grid of values of penalty parameters \eqn{\lambda}1 and \eqn{\lambda}2, stores connection signs, number of iterations until convergence and convergence outcome.
}
\usage{lasso.net.grid(x,y ,beta.0,lambda1,lambda2,M1,m.iter,n.iter,iscpp=TRUE,tol,alt.num)}
\arguments{
\item{x}{\eqn{n \times p} input data matrix}
\item{y}{response vector or size \eqn{n \times 1}}
\item{beta.0}{initial value for \eqn{\beta}. default - zero vector of size \eqn{n \times 1}}
\item{lambda1}{lasso penalty coefficient}
\item{lambda2}{network penalty coefficient}
\item{M1}{penalty matrix}
\item{m.iter}{maximum number of iterations for sign matrix updating; default - 100}
\item{n.iter}{maximum number of iterations for \eqn{\beta} updating; default - 1e5}
\item{iscpp}{binary choice for using cpp function in coordinate updates; 1 - use C++ (default), 0 - use R}
\item{tol}{convergence in \eqn{\beta} tolerance level; default - 1e-6}
\item{alt.num}{alt.num remaining iterataions are stored; default - 12}
}
\details{Fits network regression for the grid values of \eqn{\lambda}1 and \eqn{\lambda}2 using warm starts.}
\value{
\item{beta}{matrix of \eqn{\beta} coefficients, columns are for different \eqn{\lambda}1
parameters, rows \eqn{\lambda}2 parameters}
\item{mse}{mean squared error value}
\item{M}{array of connection signs. \eqn{M[,,i,j]} is the connection sign matrix for j-th \eqn{\lambda}1 value and i-th \eqn{\lambda}2 value}
\item{iterations}{matrix with stored number of steps for sign matrix to converge}
\item{update.steps}{matrix with stored number of steps for \eqn{\beta} updates to converge. (only stores the last values from connection signs iterations)}
\item{convergence.in.M}{matrix with stored values for convergence in sign matrix}
\item{convergence.in.grid}{matrix with stored values for convergence in \eqn{\beta} coefficients. If at least one \eqn{\beta} did not converge in sign matrix iterations, 0 (false) is stored, otherwise 1 (true)}
\item{xi.conv}{array with stored connection signs changes in each iteration}
\item{beta.alt}{array of coefficient vectors in case connection signs alternate}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{Weber, M., Striaukas, J., Schumacher, M., Binder, H. "Network-Constrained Covariate Coefficient and Connection Sign Estimation" (2018) <doi:10.2139/ssrn.3211163>}

\examples{
p=200
n=100
beta.0=array(1,c(p,1))
x=matrix(rnorm(n*p),n,p)
y=rnorm(n,mean=0,sd=1)
lambda1=c(0,1)
lambda2=c(0,1)
M1=diag(p)
lasso.net.grid(x, y, beta.0, lambda1, lambda2, M1)
}
