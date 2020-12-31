\name{WPP.test}
\alias{WPP.test}

\title{
 GoF tests based on the Weibull probability plot for the Weibull distribution
}
\description{
Computes the GoF tests based on the Weibull probability plot (WPP). The test statistics are similar to the coefficient of determination of the regression on the WPP: (RSB) test statistic of Smith and Bain, (REJG) test statistic of Evans, Johnson and Green, (SPP) test statistic of Coles based on the stabilized probability plot.
  
The second family includes the Shapiro-Wilk type test statistics: (SB) Shapiro and Brain and (OK) Ozturk and Korukoglu test statistic.

Two additional statistics can also be computed: the first one (ST1) is based on the kurtosis and the second one (ST2) is based on the skewness coefficient. 
}

\usage{
WPP.test(x, type = "SB", nsim = 200)
}
\arguments{
\item{x}{a numeric vector of data values.}

\item{type}{the type of the test statistic used. "OK" is the test statistic of Ozturk and Korukoglu, "RSB" the test statistic of Smith and Bain, "REJG" the test statistic of Evans-Johnson and Green, "SPP" the test statistic of Coles based on the stabilized probability plot, "SB" the test statistic of Shapiro and Brain, "ST1" the test statistic based on the kurtosis and "ST2" the test statistic based on the skewness.}

\item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}
}
\value{An object of class htest.}

\author{Meryam KRIT}
\examples{
x <- rlnorm(10,3)

#Apply the Ozturk-Korukoglu test
\dontrun{
 WPP.test(x,type="OK")
 }


# Apply the generalized smooth test based on the kurtosis 
\dontrun{
WPP.test(x,type="ST1")
}

}


\references{
\itemize{
\item Coles S.G., On goodness-of-fit tests for the two-parameter Weibull distribution derived from the stabilized probability plot, \emph{Biometrika}, 76 (3), 593-598, 1989.
  
\item Evans J.W., Johnson R.A. and Green D.W., \emph{Two and three parameter Weibull goodness-of-fit tests}, Research paper FPL-RP-493, U.S. Forest Products Laboratory, Madison, WI, 1989.
  
\item Ozturk A. and Korukoglu S., A new test for the extreme value distribution, \emph{Communications in Statistics - Simulation and Computation}, 17, 1375-1393, 1988.

\item Smith R.M. and Bain L.J., Correlation type goodness-of-fit statistics with censored sampling, \emph{Communications in statistics}, 5, 119-132, 1976.

\item Best D.J., Rayner J.C.W and Thas O., Comparison of five tests of fit for the extreme value distribution, \emph{Journal of Statistical Theory and Practice}, 1 (1), 89-99, 2007.}
}
\keyword{Stabilized Weibull probability plot}
\keyword{Shapiro-Wilk GoF tests}
\keyword{smooth tests}
\keyword{skewness}
\keyword{kurtosis}
\keyword{Weibull probability plot}
\keyword{WPP}
\keyword{Ozturk-Korukoglu test}
\keyword{Smith-Bain test}
