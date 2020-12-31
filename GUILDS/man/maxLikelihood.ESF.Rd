\name{maxLikelihood.ESF}
\alias{maxLikelihood.ESF}

\title{
Maximization of the loglikelihood given the standard Neutral Model, using the Etienne Sampling Formula
}
\description{
 This function computes the maximum likelihood estimates of the parameters of the Neutral model, using the Etienne Sampling Formula
}
\usage{
maxLikelihood.ESF(init_vals, abund, verbose = TRUE)
}

\arguments{
  \item{init_vals}{A vector of initial starting values, of the format c(theta, m)}

  \item{abund}{Vector containing a record of the number of individuals per species}

  \item{verbose}{TRUE/FALSE flag, indicates whether intermediate output is shown on screen}
}
\value{the output is a list containing the following:
	\item{par}{ a vector containing the parameter values at the maximum likelihood c(theta, m)}
	\item{fvalues}{ the likelihood at the corresponding parameter values}
	\item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence} \cr

}
\references{
	Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. Ecology Letters, 8(3), 253-260.
}

\author{
 Thijs Janzen
}

\examples{
	A <- c(1,1,1,3,5,8)
	maxLikelihood.ESF( c(7,0.1), abund = A)
}

