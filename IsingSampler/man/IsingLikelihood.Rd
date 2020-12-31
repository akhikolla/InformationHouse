\name{IsingLikelihood}
\alias{IsingLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood of all states from tractible Ising model.
}
\description{
This function returns the likelihood of all possible states. Is only tractible up to rougly 10 nodes.
}
\usage{
IsingLikelihood(graph, thresholds, beta, responses = c(0L, 1L),
                 potential = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
Square matrix indicating the weights of the network. Must be symmetrical with 0 as diagonal.
}
  \item{thresholds}{
Vector indicating the thresholds, also known as the external field.
}
  \item{beta}{
Scalar indicating the inverse temperature.
}
  \item{responses}{
Response options. Typically set to \code{c(-1L, 1L)} or \code{c(0L, 1L)} (default). Must be integers!
}
\item{potential}{
Logical, return the potential instead of the probability of each state?
}
}

\author{
Sacha Epskamp
}
