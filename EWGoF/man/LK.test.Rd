\name{LK.test}
\alias{LK.test}
\title{GoF tests based on the likelihood for the Exponential distribution

}
\description{
Computes the Exponential GoF tests based on the three following statistics: the score of Cox-Oakes (CO), Wald (W) and likelihood ratio (LR) GoF tests. These tests include the Exponential distribution in the Weibull distribution and apply a parametric test to check whether the shape parameter is equal to one.}
\usage{
LK.test(x, type = "CO", nsim = 200)
}
\arguments{
  \item{x}{a numeric vector of data values.}
  \item{type}{the type of the test statistic used. "CO" is the default used test of Cox-Oakes,"W" for Wald, "LR" for likelihood ratio test statistic.}
  \item{nsim}{an integer specifying the number of replicates used in Monte Carlo.}

}
\details{
The test statistic of Cox and Oakes is equivalent to the score test applied to the shape parameter of the Weibull distribution. The expression of the score is given using the observed information matrix and not the exact one. 
}
\value{An object of class htest.}
\references{

Cox D.R. and Oakes D., \emph{Analysis of survival data}, Chapman and Hall, 1984.

}
\author{Meryam KRIT}

\examples{
x <- rgamma(20,0.3)

#Apply the Cox and Oakes test (score) test
LK.test(x,type="CO")

# Apply the test based on Wald
LK.test(x,type="W")

# Apply the test based on the likelihood ratio
LK.test(x,type="LR")
}
\keyword{Score test}
\keyword{Wald test}
\keyword{Likelihood ratio test}
\keyword{Cox-Oakes}
