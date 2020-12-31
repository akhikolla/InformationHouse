\name{cauchypoly}
\alias{qcauchypoly}
\alias{pcauchypoly}
\alias{dcauchypoly}
\alias{rcauchypoly}
\alias{cauchypoly_inv}
\alias{cauchypoly_cdf}
\alias{cauchypoly_pdf}
\alias{cauchypoly_rnd}
\concept{quantile mixture}
\concept{distribution family}
\concept{asymmetric distribution}
\concept{heavy tails}
\title{Cauchy-polynomial quantile mixture}
\description{Density, distribution function, quantile function and random
     generation for the Cauchy-polynomial quantile mixture.}
\usage{
dcauchypoly(x,param)
pcauchypoly(x,param)
qcauchypoly(cp,param)
rcauchypoly(n,param)
cauchypoly_pdf(x,param)
cauchypoly_cdf(x,param)
cauchypoly_inv(cp,param)
cauchypoly_rnd(n,param)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{cp}{vector of probabilities}  	
  \item{n}{number of observations}
  \item{param}{vector of parameters}	
}
\value{'dcauchypoly' gives the density, 'pcauchypoly' gives the cumulative distribution
     function, 'qcauchypoly' gives the quantile function, and 'rcauchypoly'
     generates random deviates.
}
\details{The length the parameter vector specifies the order of the polynomial in the quantile mixture. 
If k<-length(param) then param[1:(k-1)] contains the mixture coefficients of polynomials starting from the constant and
param[k] is the mixture coefficient for Cauchy distribution. (Functions cauchypoly\_pdf, cauchypoly\_cdf, cauchypoly\_inv and 
cauchypoly\_rnd are aliases for compatibility with older versions of this package.)
}
\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.
}
\seealso{\code{\link{data2cauchypoly4}} for the parameter estimation and
\code{\link{dnormpoly}} for the normal-polynomial quantile mixture.}
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

