\name{IsingStateProb}
\alias{IsingStateProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood of single state from tractible Ising model.
}
\description{
This function returns the likelihood of a single possible state. Is only tractible up to rougly 10 nodes.
}
\usage{
IsingStateProb(s, graph, thresholds, beta, responses = c(0L, 1L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
Vector contaning the state to evaluate.
}
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
}

\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
