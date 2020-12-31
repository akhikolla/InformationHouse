\name{EDF_NS.test}
\alias{EDF_NS.test}
\title{
GoF tests based on the empirical distribution function, the normalized spacings and the probability plots for the Exponential distribution
}
\description{
Computes the Exponential GoF tests based on the empirical distribution function: the Kolmogorov-Smirnov (KS), Cramer-Von-Mises (CM) and Anderson-Darling (AD) tests, the tests based on the probability plot: Shapiro-Wilk (SW) and Patwardhan (PA) tests and the tests based on the normalized spacings: Gnedenko (Gn) and Gini (G) tests.
}
\usage{
EDF_NS.test(x, type = "AD", nsim = 200)
}
\arguments{
  \item{x}{a numeric vector of data values.}
  \item{type}{the type of the test statistic used. "AD" is the default used test of Anderson-Darling,"KS" for Kolmogorov-Smirnov, "CM" for Cramer-Von-Mises, "SW" for Shapiro-Wilk, "PA" for Patwardhan, "Gn" for Gnedenko and "G" for Gini test statistic.}
  \item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}
}
\details{
This function computes the GoF test statistics of three different families: the tests based on the empirical distribution function, the tests based on the probability plots and the tests based on the normalized spacings. The p-value of the tests is computed using Monte-Carlo simulations because only the asymptotic distribution of the previous statistics is known. Therefore the tests can be applied to small samples.  }

\value{An object of class htest.}

\references{
\itemize{
\item D'Agostino R.B. and Stephens M.A., \emph{Goodness-of-fit techniques}, Marcel Dekker, 1986.
\item Gail M.H. and Gastwirth J.L., A scale-free goodness-of-fit test for the exponential distribution based on the Gini statistic, \emph{Journal of the Royal Statistical Society, Series B}, 40, 350-357, 1978.
\item  Gnedenko B.V., Belyayev Y.K. and Solovyev A.D., \emph{Mathematical Models of Reliability Theory}, Academic Press, 1969.
\item Shapiro S.S. and Wilk M.B., An analysis of variance test for the exponential distribution (complete samples), \emph{Technometrics}, 14, 355-370, 1972.
\item Patwardhan G., Tests for exponentiality, \emph{Communications in Statistics, Theory and Methods}, 17, 3705-3722, 1988.
}
}
\author{
Meryam KRIT}


\examples{

x1 <- rexp(50,2)

#Apply the Kolmogorov-Smirnov test 
EDF_NS.test(x1,type="KS")

x2 <- rlnorm(50,0.2)

#Apply the Patwardhan test 
EDF_NS.test(x2,type="PA")

#Apply the Cramer-von Mises test 
EDF_NS.test(x2,type="CM")

#Apply the Gini test  
EDF_NS.test(x2,type="G")

}

\keyword{Empirical distribution function}
\keyword{Gini}
\keyword{Gndenko}
\keyword{Shapiro-Wilk}
\keyword{Patwardhan}
