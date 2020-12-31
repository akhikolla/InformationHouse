\name{normpoly}
\alias{qnormpoly}
\alias{pnormpoly}
\alias{dnormpoly}
\alias{rnormpoly}
\alias{normpoly_inv}
\alias{normpoly_cdf}
\alias{normpoly_pdf}
\alias{normpoly_rnd}
\concept{quantile mixture}
\concept{distribution family}
\concept{asymmetric distribution}
\title{Normal-polynomial quantile mixture}
\description{Density, distribution function, quantile function and random
     generation for the normal-polynomial quantile mixture.}
\usage{
dnormpoly(x,param)
pnormpoly(x,param)
qnormpoly(cp,param)
rnormpoly(n,param)
normpoly_pdf(x,param)
normpoly_cdf(x,param)
normpoly_inv(cp,param)
normpoly_rnd(n,param)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{cp}{vector of probabilities}  	
  \item{n}{number of observations}
  \item{param}{vector of parameters}	
}
\value{'dnormpoly' gives the density, 'pnormpoly' gives the cumulative distribution
     function, 'qnormpoly' gives the quantile function, and 'rnormpoly'
     generates random deviates.
}
\details{The length the parameter vector specifies the order of the polynomial in the quantile mixture. 
If k<-length(param) then param[1:(k-1)] contains the mixture coefficients of polynomials starting from the constant and
param[k] is the mixture coefficient for normal distribution. (Functions normpoly\_pdf, normpoly\_cdf, normpoly\_inv and 
normpoly\_rnd are aliases for compatibility with older versions of this package.)
}
\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.
}
\seealso{\code{\link{data2normpoly4}} for the parameter estimation and
\code{\link{dcauchypoly}} for the Cauchy-polynomial quantile mixture.
}
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
hist(x,30,freq=FALSE)
plotpoints<-seq(min(x)-1,max(x)+1,by=0.01);
lines(plotpoints,dnormpoly(plotpoints,estim_params),col='red');
lines(plotpoints,dnormpoly(plotpoints,true_params),col='blue');
}
\author{Juha Karvanen \email{juha.karvanen@iki.fi}}
\keyword{distribution}
\keyword{robust}

