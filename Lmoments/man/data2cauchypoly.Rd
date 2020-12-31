\name{data2cauchypoly}
\alias{data2cauchypoly4}
\alias{t1lmom2cauchypoly4}
\title{Estimation of the Cauchy-polynomial quantile mixture}
\description{
Estimates the parameters of the Cauchy-polynomial quantile mixture from data or from trimmed L-moments}
\usage{
data2cauchypoly4(data)
t1lmom2cauchypoly4(t1lmom)
}
\arguments{
  \item{data}{vector}
  \item{t1lmom}{vector of trimmed L-moments}
}
\value{vector containing the four parameters of the Cauchy-polynomial quantile mixture}

\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.
}
\seealso{\code{\link{t1lmoments}} for trimmed L-moments, 
\code{\link{dcauchypoly}} for the Cauchy-polynomial quantile mixture and
\code{\link{data2normpoly4}} for the estimation of the normal-polynomial quantile mixture.}
\examples{
#Generates 500 random variables from the Cauchy-polynomial quantile mixture, 
#calculates the trimmed L-moments,
#estimates parameters via trimmed L-moments and 
#plots the true pdf and the estimated pdf together with the histogram of the data.
true_params<-t1lmom2cauchypoly4(c(0,1,0.075,0.343));
x<-rcauchypoly(500,true_params);
t1lmom<-t1lmoments(x);
estim_params<-t1lmom2cauchypoly4(t1lmom);
plotpoints<-seq(-10,10,by=0.01);
histpoints<-c(seq(min(x)-1,-20,length.out=50),seq(-10,10,by=0.5),seq(20,max(x)+1,length.out=50));
hist(x,breaks=histpoints,freq=FALSE,xlim=c(-10,10));
lines(plotpoints,dcauchypoly(plotpoints,estim_params),col='red');
lines(plotpoints,dcauchypoly(plotpoints,true_params),col='blue');
}
\author{Juha Karvanen \email{juha.karvanen@iki.fi}}
\keyword{distribution}
\keyword{robust}


