\name{graphicalVARsim}
\alias{graphicalVARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates data from the graphical VAR model
}
\description{
Simulates data from the graphical VAR model, see \code{\link{graphicalVAR}} for details.
}
\usage{
graphicalVARsim(nTime, beta, kappa, mean = rep(0, ncol(kappa)), init =
                 mean, warmup = 100, lbound = rep(-Inf, ncol(kappa)),
                 ubound = rep(Inf, ncol(kappa)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nTime}{
Number of time points to sample
}
  \item{beta}{
The Beta matrix to use
}
  \item{kappa}{
The Kappa matrix to use
}
  \item{mean}{
Means to use
}
  \item{init}{
Initial values
}
  \item{warmup}{
The amount of samples to use as warmup (not returned)
}
  \item{lbound}{
Lower bound, at every time point values below this bound are set to the bound.
}
  \item{ubound}{
Upper bound, at every time point values above this bound are set to the bound.
}
}

\value{
A matrix containing the simulated data.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
