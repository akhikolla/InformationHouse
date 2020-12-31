\name{Student}
\alias{Student}
\docType{data}
\title{A simulated dataset}
\description{
  This simulated data set consists in 50 subgroups, each with 5 observations,
  on 4 variables.

  There is an isolated location shift involving
  only the first variable at time \eqn{t=10} and
  a step shift, involving the third and fourth variables,
  starting from \eqn{t=31}. The in-control distribution
  is Student's t with 3 degrees of freedom, zero mean
  and such that \eqn{cov(X_i,X_j)=0.8^{|i-j|}}{cov(Xi,Xj)=0.8^abs(i-j)}.

  See the example for the exact code used to simulate the data.
}
\usage{data(Student)}
\format{
  A 4x5x50 array.
}
\examples{
data(Student)
mphase1(Student)
#
# Replication of the simulation
#
# Generation of the in-control observations
set.seed(1)
m <- 50
n <- 5
p <- 4
df <- 3
Sigma <- outer(1:p,1:p,function(i,j) 0.8^abs(i-j))
Sigma
xnorm <- crossprod(chol(Sigma),matrix(rnorm(p*n*m),p))
xchisq <- sqrt(rchisq(n*m,df)/(df-2))
x <- array(sweep(xnorm,2,xchisq,"/"),c(p,n,m))
# Then, we add an isolated shift at time 10
# (only for the first variable)
x[1,,10] <- x[1,,10]+1
# and, a step shift starting at time 31
# (only for the third and fourth variable)
x[3:4,,31:50] <- x[3:4,,31:50] + c(0.50,-0.25)
dimnames(x)<-list(paste("X",1:4,sep=""),NULL,NULL)
identical(x,Student)
}
\keyword{datasets}
