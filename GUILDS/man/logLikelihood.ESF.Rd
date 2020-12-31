\name{logLikelihood.ESF}
\alias{logLikelihood.ESF}

\title{
	Likelihood of the Etienne sampling formula
}
\description{
	This function calculates the likelihood of the Etienne Sampling Formula, provided abundance data and parameter values.
}
\usage{
	logLikelihood.ESF(theta, m, abund)
}

\arguments{
  \item{theta}{Parameter value for the fundamental biodiversity number theta}

  \item{m}{Parameter value for migration}

  \item{abund}{Vector containing abundance data}
}
\value{Returns the LogLikelihood}

\references{
	Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. Ecology Letters, 8(3), 253-260.
}

\author{
	Thijs Janzen
}

\examples{
	A <- c(1,1,1,3,5,8); #Artificial abundance dataset
	LL <- logLikelihood.ESF(theta = 7, m = 0.1, abund = A)
}