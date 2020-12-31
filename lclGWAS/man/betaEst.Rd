\name{betaEst}
\alias{betaEst}
\docType{package}
\title{
	Estimate the Fixed Effect Parameter for Discrete-Time Multivariate Frailty Model for Grouped Survival Data
}
\description{
  A method to estimate the fixed effect parameter for a multivariate frailty model accounting for family structure of related individuals (i.e., trios). The input data is assumed to be organized such that records for each family occur consecutively, and that records for offspring precede those for parents. The variance matrix for the random effects is assumed to be of the form \code{var}*K, where K is a matrix of kinship coefficients between family members. The following groupings are permitted: (Individual), (Offspring, Offspring), (Offspring, Parent), (Offspring, Parent, Parent), and (Offspring, Offspring, Parent, Parent). Other family structures have not been implemented. 
}

\usage{
  betaEst(fam_group, alpha, dtime, delta, g, var, lower, upper)  	
}

\arguments{
  \item{fam_group}{Vector of family indicators for each sample.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{dtime}{Vector of observed discrete survival times for each sample.}
  \item{delta}{Event indicator vector: 1 indicates observed event, 0 indicates censored.}
  \item{g}{Vector of numeric genotypes for each sample.}
  \item{var}{Scalar for frailty variance.}
  \item{lower}{Scalar for the lower bound of the fixed effect parameter estimation search region.}
  \item{upper}{Scalar for the upper bound of the fixed effect parameter estimation search region.}
}


\value{
  A \code{list} with one object:
  \item{betaEst}{Scalar estimate of the fixed effect parameter.}
}


\examples{
# Generate dummy data	
fam_group <- rep(1:3,each=3)
alpha <- c(0.7500000, 0.6666667, 0.5000000, 0.0000000)
dtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)
g     <- c(0, 1, 1, 1, 2, 2, 0, 0, 0)
var   <- 0.2
res   <- betaEst(fam_group, alpha, dtime, delta, g, var, lower = 0, upper = 2) 
res
}

