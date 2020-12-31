\name{WLK.test}
\alias{WLK.test}
\title{GoF tests based on the likelihood for the Weibull distribution}
\description{
Computes the Weibull GoF tests based on the three following statistics: the score, Wald and likelihood ratio GoF tests. These tests include the Weibull distribution in larger statistics and apply a parametric test to the additional parameter.
}
\usage{
WLK.test(x, type = "GG1", funEstimate = "MLE", procedure = "S", nsim = 200, r=0)}

\arguments{
\item{x}{a numeric vector of data values.}

\item{type}{the type of the test statistic used: \itemize{
 \item "GG1" is the default used test based on the Generalized Gamma distribution, 
 \item "GG2" is also based on the Generalized Gamma distribution after a transformation,
 \item "EW" for the tests based on the Exponentiated Weibull, 
 \item "PGW" for the tests based on the Power Generalized Weibull, 
 \item "MO" for the tests based on the Marshall-Olkin distribution,
 \item "MW" for the tests based on the Modified Weibull distribution,
 \item "T1" for a combination of two tests "PGW" and "MW" (max),
 \item "T2" for a combination of two tests "PGW" and "MW" (sum),
 \item "G" for a combination of two tests based on "GG1"  (Wald and likelihood ratio).}
  }
\item{funEstimate}{the method used to estimate the two Weibull parameters. "MLE" is the default used method based on the Maximum Likelihood Estimators, "LSE" for the Least Squares Estimators and "ME" for the Moment Estimators.}

\item{procedure}{the procedure used as a default is the score "S". The procedure can be either "W" for the Wald test or "LR" for the test based on the likelihood ratio procedure.}

\item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}

\item{r}{an integer specifying the number of right censored observations.}
}
\value{An object of class htest.}
\details{

The tests are based on different generalized Weibull families: the tests (GG1) and (GG2) are based on the Generalized Gamma distribution, the tests (EW) are based on the Exponentiated Weibull, (PGW) on the Power Generalized Weibull, (MO) on the Marshall-Olkin distribution and (MW) are based on the Modified Weibull distribution. Each family can have nine versions depending on the procedure used (score, Wald or likelihood ratio statistic) and on the parameters estimation methods: maximum likelihood, moment or least squares method, except GG1 which has only three versions using the maximum likelihood estimators. 

The tests statistics \eqn{T1} and \eqn{T2} are a combination between two Wald test statistics: PGW with ME (\eqn{\breve{PGW}_w}) and MW with MLE (\eqn{\hat{MW}_w}) after they are centered with their mean values (\eqn{\overline{\breve{PGW}}_w} and \eqn{\overline{\hat{MW}}_w}) and normalized by their standard deviations (respectively \eqn{sd(\breve{PGW}_w)} and \eqn{sd(\hat{MW}_w)}). 

The expressions of the statistics \eqn{T1} and \eqn{T2} are as follows:

\eqn{T1=max(\left|\frac{\breve{PGW}_w-\overline{\breve{PGW}_w}}{sd(\breve{PGW}_w)}\right|,\left|\frac{\hat{MW}_w-\overline{\hat{MW}_w}}{sd(\hat{MW}_w)}\right|)}

\eqn{T2=0.5\left|\frac{\breve{PGW}_w-\overline{\breve{PGW}_w}}{sd(\breve{PGW}_w)}\right|+0.5\left|\frac{\hat{MW}_w-\overline{\hat{MW}_w}}{sd(\hat{MW}_w)}\right|}

All the previous tests can be applied to type II right censored samples (simple censoring). The censoring is introduced in the MLEs. 
A second statistic \eqn{G} similar to \eqn{T}, is combining Wald and likelihood ratio tests based on the Generalized Gamma distribution:

\eqn{G=0.5\left|\frac{\hat{GG}^1_l-\overline{\hat{GG}^1_l}}{sd(\hat{GG}^1_l)}\right|+0.5\left|\frac{\hat{GG}^1_w-\overline{\hat{GG}^1_w}}{sd(\hat{GG}^1_w)}\right|}

}
\references{
Krit M., Gaudoin O., Xie M. and Remy E., Simplified likelihood goodness-of-fit tests for the Weibull distribution, \emph{Communications in Statistics - Simulation and Computation}.}
\author{Meryam KRIT}
\examples{

x <- rlnorm(50,.3)

#Apply some likelihood based tests
WLK.test(x,type="GG1",funEstimate="MLE",procedure="W")
WLK.test(x,type="PGW",funEstimate="ME",procedure="S")
WLK.test(x,type="MO",funEstimate="LSE",procedure="LR")

#Apply G to censored sample at right r=10
a<- sort(x[1:40])
\dontrun{
WLK.test(a,type="G",r=10)
}
}

\keyword{Likelihood based tests}
\keyword{Score statistic}
\keyword{Wald statistic}
\keyword{Likelihood ratio statistic}

