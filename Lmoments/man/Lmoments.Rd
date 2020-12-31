\name{Lmoments}
\alias{Lmoments}
\alias{Lcoefs}
\alias{Lmomcov}
\alias{Lmoments_calc}
\alias{Lmomcov_calc}
\alias{shiftedlegendre}
\concept{moments}
\concept{skewness}
\concept{kurtosis}
\concept{descriptive statistics}
\title{L-moments}
\description{
  Calculates sample L-moments, L-coefficients and covariance matrix of L-moments. 
}
\usage{
Lmoments(data, rmax = 4, na.rm = FALSE, returnobject = FALSE, trim = c(0, 0))
Lcoefs(data, rmax = 4, na.rm = FALSE, trim = c(0, 0))
Lmomcov(data, rmax = 4, na.rm = FALSE)
Lmoments_calc(data, rmax = 4)
Lmomcov_calc(data, rmax = 4)
shiftedlegendre(rmax)
}
\arguments{
  \item{data}{matrix or data frame.}
  \item{rmax}{maximum order of L-moments.}
  \item{na.rm}{a logical value indicating whether 'NA' values should be
          removed before the computation proceeds.}
  \item{returnobject}{a logical value indicating whether a list object should be returned instead of an array of L-moments.}
  \item{trim}{c(0, 0) for ordinary L-moments and c(1, 1) for trimmed (t = 1) L-moments}
}
\value{
\code{Lmoments} returns an array of L-moments containing a row for each variable in data, or if returnobject=TRUE, 
a list containing 
\item{lambdas}{an array of L-moments} 
\item{ratios}{an array of mean, L-scale and L-moment ratios}
\item{trim}{the value of the parameter 'trim'}
\item{source}{a string with value "Lmoments" or "t1lmoments".}

\code{Lcoefs} returns an array of L-coefficients (mean, L-scale, L-skewness, L-kurtosis, ...) 
containing a row for each variable in data.

\code{Lmomcov} returns the covariance matrix of L-moments or a list of covariance matrices if the input has multiple columns.
The numerical accuracy of the results decreases with increasing \code{rmax}.
With \code{rmax > 5}, a warning is thrown, as the numerical accuracy of the results is likely less than \code{sqrt(.Machine$double.eps)}. 

\code{shiftedlegendre} returns a matrix of the coefficients of the shifted Legendre polynomials up to a given order.}

\note{Functions \code{Lmoments} and \code{Lcoefs} calculate trimmed L-moments if you specify \code{trim = c(1, 1)}.
\code{Lmoments_calc} and \code{Lmomcov_calc} are internal C++ functions called by \code{Lmoments} and \code{Lmomcov}.
The direct use of these functions is not recommended.}

\references{Karvanen, J. 2006. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
\emph{Computational Statistics & Data Analysis} \bold{51}, (2), 947--959.
\url{http://www.bsp.brain.riken.jp/publications/2006/karvanen_quantile_mixtures.pdf}.

Elamir, E. A., Seheult, A. H. 2004. Exact variance structure of sample L-moments, \emph{Journal of 
Statistical Planning and Inference} \bold{124} (2) 337--359.

Hosking, J. 1990. L-moments: Analysis and estimation distributions using linear combinations of order statistics, 
\emph{Journal of Royal Statistical Society B} \bold{52}, 105--124.
}
\seealso{\code{\link{t1lmoments}} for trimmed L-moments,
\code{\link{dnormpoly}}, \code{\link{lmom2normpoly4}} and \code{\link{covnormpoly4}} for the normal-polynomial quantile mixture
and package \pkg{lmomco} for additional L-moment functions}
\examples{
#Generates a sample 500 observations from the normal-polynomial quantile mixture, 
#calculates the L-moments and their covariance matrix,
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
\author{Juha Karvanen \email{juha.karvanen@iki.fi}, Santeri Karppinen}
\keyword{robust}
\keyword{univar}

