\name{t1lmoments}
\alias{t1lmoments}
\alias{t1lmoments_calc}
\concept{moments}
\concept{skewness}
\concept{kurtosis}
\concept{descriptive statistics}
\concept{heavy tails}
\title{Trimmed L-moments}
\description{
Calculates sample trimmed L-moments with trimming parameter 1. }
\usage{
t1lmoments(data, rmax = 4)
t1lmoments_calc(data, rmax = 4)
}
\arguments{
  \item{data}{matrix or data frame.}
    \item{rmax}{maximum order of trimmed L-moments.}
}
\value{array of trimmed L-moments (trimming parameter = 1) up to order 4 containing a row for each variable in data.}

\note{Functions \code{link{Lmoments}} and \code{link{Lcoefs}} calculate trimmed L-moments if you specify \code{trim = c(1, 1)}.
\code{t1lmoments_calc} is an internal C++ function called by \code{t1lmoments}. The direct use of this function is not recommended.}

\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.

Elamir, E. A., Seheult, A. H. 2003. Trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{43}, 299--314.
}
\seealso{\code{\link{Lmoments}} for L-moments, and 
\code{\link{dcauchypoly}} and \code{\link{t1lmom2cauchypoly4}} for the Cauchy-polynomial quantile mixture}
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
\author{Juha Karvanen \email{juha.karvanen@iki.fi}, Santeri Karppinen}
\keyword{robust}
\keyword{univar}


