\name{fastols}
\alias{fastols}
\title{
Fast least squares estimate.

}
\description{
Computes least squares estimate in an efficient way.
}
\usage{fastols(y, x)}
\arguments{
\item{y}{dependent variable}
\item{x}{response variable}
}
\author{
Maintainer: Jonas Striaukas <jonas.striaukas@gmail.com>
}


\examples{
p<-10
n<-100
x<-matrix(rnorm(n*p),n,p)
beta<-array(5, c(p,1))
y<-x\%*\%beta + rnorm(n,mean=0,sd=0.1)
fastols(y,x)
}
