\name{WLP.test}
\alias{WLP.test}
\title{
GoF tests based on the Laplace transform for the Weibull distribution
}
\description{
Computes the Weibull GoF tests based on the Laplace transform: (CQ) Cabana-Quiroz with the two versions either with the Maximum Likelihood Estimators (MLE) or the Moment Estimators (ME) and (LT) based on the Laplace transform with discretization which has three versions, it can be computed either with the MLE, ME or the Least Squares Estimators (LSE). 
}
\usage{WLP.test(x,type = "CQ",funEstimate="MLE",s1=-1,s2=0.4,mr=100,a=-5,nsim = 200)}
\arguments{
\item{x}{a numeric vector of data values.}

\item{type}{the type of the test statistic used. "CQ" is the default used test of Cabana and Quiroz and "LT" is the test statistic based on the Laplace transform discretization.}

\item{funEstimate}{the method used to estimate the two Weibull parameters. "MLE" is the default used method based on the maximum likelihood estimators, "LSE" for the least squares estimators and "ME" for the moment estimators. The test statistic "CQ" can be computed for MLE or ME only.} 

\item{s1}{the first value where the quadratic form of the statistic "CQ" is computed.}

\item{s2}{the second value where the quadratic form of the statistic "CQ" is computed.}
  
\item{mr}{the number of the discretizations done to compute the test "LT"; mr = 100 is the default value recommended.}

\item{a}{the value of the adjusting parameter used in the statistic "LT"; a = -5 is the default value used.}
    
\item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}
}
\value{An object of class htest.}
\references{
Cabana A. and Quiroz A.J., Using the empirical moment generating function in testing the Weibull and type 1 Extreme Value distributions, \emph{ Test}, 14(2), 417-431, 2005.

Krit M., Goodness-of-fit tests for the Weibull distribution based on the Laplace transform, \emph{J-SFdS}, 2014.}

\author{Meryam KRIT}
\examples{
# Xie data
x<-c(5, 11, 21, 31, 46, 75, 98, 122, 145, 165, 195, 224, 245, 293, 321, 330, 350, 420)

#Apply test of Cabana and Quiroz with the least squares estimators
\dontrun{WLP.test(x, type="CQ*",funEstimate="LSE")}
}
\keyword{Laplace transform}
\keyword{Cabana-Quiroz test}
