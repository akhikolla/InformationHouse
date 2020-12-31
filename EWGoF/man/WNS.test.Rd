\name{WNS.test}
\alias{WNS.test}
\title{GoF tests based on the normalized spacings for the Weibull distribution}
\description{
Computes the Weibull GoF tests based on the normalized spacings: (TS) Tiku-Singh, (MSF) Mann-Scheuer-Fertig only with right censoring and (LOS) Lockhart-O'Reilly-Stephens test statistics.
}
\usage{
WNS.test(x, type = "TS", s = 0, r = 0, nsim = 200)
}
\arguments{
\item{x}{a numeric vector of data values.}

\item{type}{the type of the test statistic used. "TS" is the default used test of Tiku-Singh,"MSF" for Mann-Scheuer-Fertig and "LOS" for Lockhart-O'Reilly-Stephens test statistic.}

\item{s}{the index of the smallest observed value of x (s=0 if no left censoring).}

\item{r}{the index of the largest observed value of x (r=0 if no right censoring).}

\item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}
}
\details{
For these tests statistics, there is no need to estimate the unknown values of the Weibull parameters and they are among the most powerful ones especially TS and LOS.

The p-value computed is not the exact p-value: the null hypothesis distribution of some statistics is known only asymptotically and sometimes is not even known. The asymptotic approximation is not correct especially for small samples. That is why Monte Carlo simulation is needed to compute the p-value.
}
\value{An object of class htest.}
\references{
Tiku M.L. and Singh M., Testing the two-parameter Weibull distribution, \emph{Communications in Statistics}, 10, 907-918, 1981.

Mann N.R., Scheuer E.M. and Fertig K.W., A new goodness-of-fit test for the two-parameter Weibull or extreme-value distribution, \emph{Communications in Statistics}, 2, 383-400, 1973.

Lockhart R.A., O'Reilly F. and Stephens M.A., Tests for the extreme-value and Weibull distributions based on normalized spacings, \emph{Naval Research Logistics Quarterly}, 33, 413-421, 1986.
}
\author{Meryam KRIT
}
\examples{
x1 <- rweibull(50,2,3)

#Apply the Tiku-Singh test with censored sample: 4 left censored and 5 right censored values
\dontrun{
WNS.test(x1,type="TS",s=4,r=5)
}



}
\keyword{normalized spacings}
\keyword{Tiku-Singh statistic}
\keyword{Mann-Scheuer-Fertig statistic}
\keyword{Lockhart-O'Reilly-Stephens test statistic}
