\name{alphaEst}
\alias{alphaEst}
\docType{package}
\title{
	Estimate the Interval Baseline Survival Rates for Discrete-Time Multivariate Frailty Model for Grouped Survival Data
}
\description{
  A method to estimate the baseline survival rate for each time interval for a multivariate frailty model. The estimation is conduced under the null hypothesis, i.e., that there is no fixed effect.
}

\usage{
  alphaEst(dtime, delta)  	
}

\arguments{
  \item{dtime}{Vector of observed discrete survival times for each sample.}
  \item{delta}{Event indicator vector: 1 indicates observed event, 0 indicates censored.}
}


\value{
  A \code{list} with one object:
  \item{alphaEst}{Vector of estimates of the baseline survival rates.}
}


\examples{
# Generate dummy data	
dtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)
res   <- alphaEst(dtime, delta) 
res
}

