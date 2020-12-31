\name{data2normpoly}
\alias{data2normpoly4}
\alias{lmom2normpoly4}
\alias{data2normpoly6}
\alias{lmom2normpoly6}
\title{Estimation of normal-polynomial quantile mixture}
\description{
Estimates the parameters of normal-polynomial quantile mixture from data or from L-moments}
\usage{
data2normpoly4(data)
lmom2normpoly4(lmom)
data2normpoly6(data)
lmom2normpoly6(lmom)
}
\arguments{
  \item{data}{matrix or data frame}
  \item{lmom}{vector or matrix of L-moments}
}
\value{vector or matrix containing the four or six parameters of normal-polynomial quantile mixture}
\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.
}
\seealso{\code{\link{dnormpoly}} for L-moments, 
\code{\link{dnormpoly}} for the normal-polynomial quantile mixture and
\code{\link{data2cauchypoly4}} for the estimation of Cauchy-polynomial quantile mixture.}
\examples{
#Generates a sample 500 observations from the normal-polynomial quantile mixture, 
#calculates L-moments and their covariance matrix,
#estimates parameters via L-moments and 
#plots the true pdf and the estimated pdf together with the histogram of the data.
true_params<-lmom2normpoly4(c(0,1,0.2,0.05));
x<-rnormpoly(500,true_params);
lmoments<-Lmoments(x);
lmomcov<-Lmomcov(x);
estim_params<-lmom2normpoly4(lmoments);
hist(x,30,freq=FALSE);
plotpoints<-seq(min(x)-1,max(x)+1,by=0.01);
lines(plotpoints,dnormpoly(plotpoints,estim_params),col='red');
lines(plotpoints,dnormpoly(plotpoints,true_params),col='blue');
}
\author{Juha Karvanen \email{juha.karvanen@iki.fi}}
\keyword{distribution}
\keyword{robust}

