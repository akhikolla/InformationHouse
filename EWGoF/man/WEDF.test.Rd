\name{WEDF.test}
\alias{WEDF.test}

\title{GoF tests based on the empirical distribution function for the Weibull distribution
}
\description{
Computes the Weibull GoF tests based on the empirical distribution function: Kolmogorov-Smirnov (KS), Cramer-Von-Mises (CM), Watson (W), Anderson-Darling (AD) and Liao-Shimokawa (LS) statistics using the three following estimation methods: Maximum Likelihood Estimators (MLE), Least Squares Estimators (LSE) and Moments Estimators (ME). A test statistic (KL) is added to this family but it is based on the Kullback-Leibler information; this test can also be computed with the three estimation methods.
}
\usage{
WEDF.test(x, type = "AD", funEstimate = "MLE", paramKL = 2, nsim = 200)
}

\arguments{
\item{x}{a numeric vector of data values.}
  
\item{type}{the type of the test statistic used, "AD" is the default used test of Anderson-Darling,"KS" for Kolmogorov-Smirnov, "CM" for Cramer-Von-Mises, "W" for Watson, "LS" for Liao-Shimokawa and "KL" for Kullback-Leibler.}
  
\item{funEstimate}{the method used to estimate the two Weibull parameters. "MLE" is the default used method based on the maximum likelihood estimators, "LSE" for the least squares estimators and "ME" for the moment estimators.}

\item{paramKL}{the value of the parameter m used in the expression of the statistic KL.}

\item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}
}
\details{
The above test statistics are known in the literature by using the maximum likelihood estimators except the test "LS" that Liao and Shimokawa recommend to use with the least squares estimators. Each of the tests can have three versions, depending on the estimation method used. 

All these tests statistics can be used for small samples. The asymptotic distributions of the tests are known in some cases but the use of their quantiles for small samples may lead to wrong conclusions. That is why we use Monte-Carlo simulations to apply the tests. 
}
\value{An object of class htest.}
\references{
\itemize{
\item D'Agostino R.B. and Stephens M.A., \emph{Goodness-of-fit techniques}, Marcel Dekker, 1986.

\item Liao M. and Shimokawa T., A new goodness-of-fit test for type-I extreme-value and 2-parameter Weibull distributions with estimated parameters, \emph{Journal of Statistical Computation and Simulation}, 64 (1), 23-48, 1999.

\item Chandra M., Singpurwalla N.D and Stephens M.A, Kolmogorov statistics for tests of fit for the Extreme Value and Weibull  distributions, \emph{ Journal of American Statistical Association}, 76 (375), 729-731, 1981. 

\item Vaquera-Huerta P.P. and Villasenor-Alva J.A., A Goodness-of-fit test for the Gumbel distribution based on Kullback-Leibler information, \emph{Communications in Statistics, Theory and Methods}, 38: 842-855, 2009.}
}
\author{
Meryam KRIT
}
\examples{

x1 <- rweibull(50,2,3)

#Apply Kolmogorov-Smirnov test with maximum likelihood estimators
WEDF.test(x1,type="KS",funEstimate="MLE")

x2 <- rlnorm(50,0.2)

#Apply the Liao and Shimokawa test using the least squares estimators
WEDF.test(x2,type="LS",funEstimate="LSE")

#Apply the Cramer-von Mises test with the moment estimators
WEDF.test(x2,type="CM",funEstimate="ME")

#Apply the test based on the Kullback-Leibler information with the moment estimators
WEDF.test(x2,type="KL",funEstimate="ME")
}

\keyword{Cramer-von-Mises test}
\keyword{Anderson-Darling test}
\keyword{Kolmogorov-Smirnov test}
\keyword{Liao and Shimokawa test}
\keyword{Kullback-Lebler information}
\keyword{Empirical distribution function}
