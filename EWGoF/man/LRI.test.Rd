\name{LRI.test}
\alias{LRI.test}

\title{
GoF tests based on the Laplace transform, the mean residual life and the integrated distribution function for the Exponential distribution 
}
\description{
Computes the Weibull GoF tests based on the Laplace transform: Baringhaus-Henze (BH) and Henze (He). The test statistic of Klar (Kl) is based on the integrated distribution function. Two tests are based on the mean residual life (BHC, BHK).}
\usage{
LRI.test(x, type = "BH", a = 1, nsim = 200)
}
\arguments{
 \item{x}{a numeric vector of data values.}
  \item{type}{the type of the test statistic used. "BH" is the default used test of Baringhaus-Henze,"He" for Henze, "Kl" for Klar, "BHC" and "BHK" for the tests based on the integrated distribution function.}
    \item{a}{parameter value to be adjusted for the test statistics.}
  \item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}

}

\value{An object of class htest.}
\references{
Baringhaus L. and Henze N., Tests of fit for exponentiality based on a characterization via the mean residual life function, \emph{Statistical Papers}, 41, 225-236, 2000.

Baringhaus L. and Henze N., A class of consistent tests for exponentiality based on the empirical Laplace transform, \emph{Annals of the Institute of Statistical Mathematics}, 43, 551-564, 1991.

Henze N., A new flexible class of omnibus tests for exponentiality, \emph{Communications in Statistics, Theory and Methods}, 22, 115-133, 1993.

Klar B., Goodness-of-fit tests for the exponential and normal distribution based on the integrated distribution function, \emph{Annals of the Institute of Statistical Mathematics}, 53, 338-353, 2001.

}
\author{Meryam KRIT}
\examples{
x <- rlnorm(50,0.3)

#Apply the Baringhaus-Henze test
LRI.test(x,type="BH")

# Apply the test of Henze
LRI.test(x,type="He")


# Apply the test of Klar
LRI.test(x,type="Kl")


# Apply the test of Barighaus based on the integrated distribution function
LRI.test(x,type="BHC")

}
\keyword{Laplace transform}
\keyword{integrated distribution function}
\keyword{mean residual life}
\keyword{Klar}
\keyword{Henze}
