\name{beta.update.net}
\alias{beta.update.net}
\title{
Updates \eqn{\beta} coefficients.
}
\description{
This function updates \eqn{\beta} for given penalty parameters.
}
\usage{beta.update.net(x,y,beta,lambda1,lambda2,M1,n.iter,iscpp,tol)}
\arguments{
\item{x}{input data matrix of size \eqn{n \times p}; n - number of observations; p - number of covariates}
\item{y}{response vector or size \eqn{n \times 1}}
\item{beta}{initial value for \eqn{\beta}; default - zero vector of size \eqn{n \times 1}}
\item{lambda1}{lasso penalty parameter}
\item{lambda2}{network penalty parameter}
\item{M1}{penalty matrix}
\item{n.iter}{maximum number of iterations for \eqn{\beta} step; default - 1e5}
\item{iscpp}{binary choice for using cpp function in coordinate updates; 1 - use C++ (default), 0 - use R}
\item{tol}{convergence tolerance level; default - 1e-6}
}
\details{Updates the coefficient vector \eqn{\beta} given the data and penalty parameters \eqn{\lambda}1 and \eqn{\lambda}2.
Convergence criterion is defined as \eqn{\sum_{i=1}^p |\beta_{i,j} - \beta_{i,j-1}| \leq} to.}
\value{
\item{beta}{updated \eqn{\beta} vector}
\item{convergence}{binary variable; 1 - yes}
\item{steps}{number of steps until convergence}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}
\references{Weber, M., Striaukas, J., Schumacher, M., Binder, H. "Network-Constrained Covariate Coefficient and Connection Sign Estimation" (2018) <doi:10.2139/ssrn.3211163>}

\examples{
p<-200
n<-100
beta.0=array(1,c(p,1))
x<-matrix(rnorm(n*p),n,p)
y<-rnorm(n,mean=0,sd=1)
lambda1<-1
lambda2<-1
M1<-diag(p)
updates<-beta.update.net(x, y, beta.0, lambda1, lambda2, M1)
}
