\name{IsingEntrophy}
\alias{IsingEntrophy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Entropy of the Ising Model
}
\description{
Returns (marginal/conditional) Shannon information of the Ising model.
}
\usage{
IsingEntrophy(graph, thresholds, beta = 1, conditional = numeric(0), 
          marginalize = numeric(0), base = 2, responses = c(0L, 1L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
Weights matrix
}
  \item{thresholds}{
Thresholds vector
}
  \item{beta}{
Inverse temperature
}
  \item{conditional}{
Indices of nodes to condition on
}
  \item{marginalize}{
Indices of nodes to marginalize over
}
  \item{base}{
Base of the logarithm
}
  \item{responses}{
Vector of outcome responses.
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
