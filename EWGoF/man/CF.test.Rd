\name{CF.test}
\alias{CF.test}
\title{
GoF tests based on the empirical characteristic function for the Exponential distribution}
\description{
Computes the GoF tests based on the characteristic function of the Exponential distribution: Epps-Pulley (EP), Henze-Meintanis (W1, W2) and Meintanis-Iliopoulos test statistics (T1, T2).}
\usage{
CF.test(x, type = "EP", a = 1, nsim = 200)
}
\arguments{
  \item{x}{a numeric vector of data values.}
  \item{type}{the type of the test statistic used. "EP" is the default used test of Epps-Pulley,"W1" and "W2" for Henze and Meintanis, "T1" and "T2" for Meintanis-Iliopoulos test statistics.}
    \item{a}{parameter value to be adjusted for the test statistics ("W1", "W2", "T1" and "T2").}
  \item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}

}
\details{
The computation time of this function is quite long for the test statistics "W1", "W2", "T1" and "T2" because of their complex expression. The Monte-Carlo simulations take more time compared to the other tests. These tests are not defined for a=0.}
\value{An object of class htest.}
\references{
 Epps T.W. and Pulley L.B., A test for exponentiality vs. monotone hazard alternatives derived from the empirical characteristic function, \emph{Journal of the Royal Statistical Society, Series B}, 48, 206-213, 1986.

Henze N. and Meintanis S.G., Recent and classical tests for exponentiality: partial review with comparisons, \emph{Metrika}, 61, 29-45, 2005.

Henze N. and Meintanis S.G., Goodness-of-fit tests based on a new characterization of the exponential distribution, \emph{Communications in Statistics, Theory and Methods}, 31, 1479-1497, 2002.

Meintanis S.G. and Iliopoulos G., Characterizations of the exponential distribution based on certain properties of its characteristic function,
\emph{Kybernetika}, 39 (3), 295-298, 2003.

}
\author{
Meryam KRIT}
\examples{
x <- rgamma(10,0.3)

#Apply the Epps-Pulley test
CF.test(x,type="EP")

# Apply the test of Meintanis-Iliopoulos
CF.test(x,type="T1",nsim=200)


# Apply the test of Henze-Meintanis
CF.test(x,type="W1",nsim=200)

}
\keyword{Epps-Pulley}
\keyword{Henze-Meintanis}
\keyword{Meintanis-Iliopoulos}
\keyword{Empirical characteristic function}
